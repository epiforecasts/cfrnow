# Simulate replicate outcomes over posterior draws and replay the real-time
# observation process, then aggregate. Kept free of the fitted model so it can
# be tested directly: `cfr`, `loc` and `sc` (and their recovery counterparts)
# are ndraws-by-n matrices of per-draw, per-case parameters on the response
# scale, `h` is the per-case follow-up horizon, and `observed` holds the
# observed death and recovery counts and the observed death delays.
.cfr_ppc_stats <- function(cfr, loc, sc, h, fam, use_rec, rloc, rsc, rfam,
                           observed) {
  nd <- nrow(cfr)
  n <- ncol(cfr)
  # Each family's draws come from its brms mu-form: mu is the delay's mean and
  # the second parameter is the family's own shape (sdlog for a lognormal).
  draw_delay <- function(loc_i, sc_i, family) {
    switch(family,
      lognormal = stats::rlnorm(length(loc_i), loc_i, sc_i),
      gamma = stats::rgamma(length(loc_i), shape = sc_i, rate = sc_i / loc_i),
      weibull = stats::rweibull(
        length(loc_i), sc_i, loc_i / gamma(1 + 1 / sc_i)
      ),
      stop("pp_check_cfr() supports lognormal, gamma and weibull delays only.",
        call. = FALSE
      )
    )
  }

  counts <- vector("list", nd)
  delays <- vector("list", nd)
  for (k in seq_len(nd)) {
    fatal <- stats::runif(n) < cfr[k, ]
    frac <- stats::runif(n) # true onset is uniform within its recorded day
    delay_day <- floor(frac + draw_delay(loc[k, ], sc[k, ], fam))
    obs_death <- fatal & (delay_day <= h - 1)

    cts <- data.frame(
      .draw = k, outcome = "deaths", n = sum(obs_death),
      stringsAsFactors = FALSE
    )
    if (use_rec) {
      rec_day <- floor(frac + draw_delay(rloc[k, ], rsc[k, ], rfam))
      obs_rec <- !fatal & (rec_day <= h - 1)
      cts <- rbind(cts, data.frame(
        .draw = k, outcome = "recoveries", n = sum(obs_rec),
        stringsAsFactors = FALSE
      ))
    }
    counts[[k]] <- cts
    dd <- delay_day[obs_death] # may be empty when a draw has no observed deaths
    delays[[k]] <- data.frame(.draw = rep(k, length(dd)), delay = dd)
  }

  observed_counts <- data.frame(
    outcome = "deaths", n = observed$deaths,
    stringsAsFactors = FALSE
  )
  if (use_rec) {
    observed_counts <- rbind(observed_counts, data.frame(
      outcome = "recoveries", n = observed$recoveries,
      stringsAsFactors = FALSE
    ))
  }

  list(
    counts = do.call(rbind, counts),
    observed_counts = observed_counts,
    delays = do.call(rbind, delays),
    observed_delays = data.frame(delay = observed$death_delays)
  )
}

# Pull the per-draw, per-case CFR and delay parameters out of the fit with
# posterior_linpred (so covariate and time-varying fits work the same as
# intercept-only ones), build the follow-up horizon, then hand off to
# .cfr_ppc_stats.
.cfr_replicate <- function(object, ndraws) {
  onset <- object$cfrnow$onset
  if (is.null(onset)) {
    stop("pp_check_cfr() needs per-case onset dates; refit from ",
      "prepare_cfr_data() output.",
      call. = FALSE
    )
  }
  d <- object$data
  if (length(onset) != nrow(d)) {
    stop("pp_check_cfr(): stored onset does not line up with the fitted data.",
      call. = FALSE
    )
  }
  fam <- object$cfrnow$family
  scale_dpar <- if (fam == "lognormal") "sigma" else "shape"

  total <- posterior::ndraws(object)
  ids <- if (ndraws >= total) {
    seq_len(total)
  } else {
    sort(sample.int(total, ndraws))
  }
  lp <- function(dpar) {
    brms::posterior_linpred(
      object,
      dpar = dpar, transform = TRUE, draw_ids = ids
    )
  }
  cfr <- lp("prob")
  loc <- lp("mu")
  sc <- lp(scale_dpar)

  use_rec <- isTRUE(object$cfrnow$use_recovery)
  rloc <- rsc <- rfam <- NULL
  if (use_rec) {
    rfam <- object$cfrnow$recovery_family
    rscale_dpar <- if (rfam == "lognormal") "rsigma" else "rshape"
    rloc <- lp("rmu")
    rsc <- lp(rscale_dpar)
  }

  # Per-case follow-up horizon: days from the recorded onset to the cut-off. NA
  # obs_time (a retrospective fit) means no truncation.
  obs_time <- object$cfrnow$obs_time %||% as.Date(NA)
  h <- if (is.na(obs_time)) {
    rep(Inf, nrow(d))
  } else {
    as.numeric(as.Date(obs_time) - as.Date(onset)) + 1
  }

  observed <- list(
    deaths = sum(d$outcome == .CURE_DEATH),
    recoveries = if (use_rec) sum(d$outcome == .CURE_RECOVERY) else NA_integer_,
    death_delays = d$y[d$outcome == .CURE_DEATH]
  )
  .cfr_ppc_stats(cfr, loc, sc, h, fam, use_rec, rloc, rsc, rfam, observed)
}

.cfr_ppc_counts_plot <- function(reps) {
  ggplot2::ggplot(reps$counts, ggplot2::aes(x = .data[["n"]])) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.6) +
    ggplot2::geom_vline(
      data = reps$observed_counts,
      ggplot2::aes(xintercept = .data[["n"]]), linewidth = 1
    ) +
    ggplot2::facet_wrap(~outcome, scales = "free") +
    ggplot2::labs(
      x = "count observed by the cut-off",
      y = "posterior-predictive draws",
      title = "Posterior-predictive check: observed counts",
      subtitle = "histogram = replicates, line = observed"
    ) +
    ggplot2::theme_minimal()
}

.cfr_ppc_delay_plot <- function(reps) {
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["delay"]])) +
    ggplot2::stat_ecdf(
      data = reps$delays, ggplot2::aes(group = .data[[".draw"]]),
      geom = "step", colour = "steelblue", alpha = 0.2
    ) +
    ggplot2::stat_ecdf(
      data = reps$observed_delays, geom = "step", linewidth = 1
    ) +
    ggplot2::labs(
      x = "observed onset-to-death delay (days)", y = "ECDF",
      title = "Posterior-predictive check: observed delays",
      subtitle = "blue = replicates, black = observed"
    ) +
    ggplot2::theme_minimal()
}

#' Posterior-predictive check for a cfrnow fit
#'
#' Draws replicate line-list outcomes from the posterior and compares them with
#' the data the model was fit to, replaying the real-time observation process
#' (the onset-to-event timing and the cut-off censoring). Two checks are
#' available: the count of observed deaths by the cut-off (plus recoveries for a
#' two-outcome fit), and the distribution of the observed onset-to-death delays.
#'
#' The check reuses the fit's own posterior draws of the CFR and the delay, so
#' it works for covariate and time-varying `prob ~ ...` fits as well as
#' intercept-only ones. It needs the observation cut-off, which [fit_cfr()]
#' records when the data come from [prepare_cfr_data()]; a retrospective fit
#' (`obs_time = NULL`) has no truncation to replay, so every fatal case shows up
#' as a death.
#'
#' Only the quantities the model generates are checked: the observed death (and
#' recovery) counts and the observed onset-to-death delays. The split of the
#' remaining cases into censored versus untimed-resolved is not part of the
#' generative model, so it is left out.
#'
#' The check is stochastic: it subsamples posterior draws and simulates
#' outcomes, so call [set.seed()] beforehand for reproducible plots.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @param type Which check to plot: `"counts"` (default) or `"delay"`.
#' @param ndraws Number of posterior draws to replicate over. Defaults to 100,
#'   or fewer if the fit has fewer draws.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' fit <- fit_cfr(d, delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)))
#' pp_check_cfr(fit, type = "counts")
#' pp_check_cfr(fit, type = "delay")
#' }
#' @family fit
#' @export
pp_check_cfr <- function(object, type = c("counts", "delay"), ndraws = 100) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("pp_check_cfr() needs the ggplot2 package.", call. = FALSE)
  }
  type <- match.arg(type)
  reps <- .cfr_replicate(object, ndraws)
  if (type == "counts") {
    .cfr_ppc_counts_plot(reps)
  } else {
    .cfr_ppc_delay_plot(reps)
  }
}
