#' Naive deaths / cases ratio
#'
#' @param n_deaths,n_cases Death and case counts.
#' @return The naive ratio, or `NA` when there are no cases.
#' @noRd
naive_cfr <- function(n_deaths, n_cases) {
  if (n_cases == 0) NA_real_ else n_deaths / n_cases
}

#' Mean and sd (days) of a delay from its brms native parameters
#'
#' `loc_var` is the location intercept and `scale_var` the second parameter's
#' intercept; the second parameter is log-linked. Handles lognormal (loc is
#' meanlog), gamma and weibull (loc is log-mean, second parameter is the shape).
#' @noRd
.delay_moments <- function(dr, loc_var, scale_var, family) {
  loc <- dr[[loc_var]]
  sc <- exp(dr[[scale_var]])
  if (family == "lognormal") {
    dmean <- exp(loc + sc^2 / 2)
    list(mean = dmean, sd = sqrt(exp(sc^2) - 1) * dmean)
  } else if (family == "weibull") {
    # brms weibull: mu is the mean (log link), sc is the shape k
    dmean <- exp(loc)
    cv <- sqrt(gamma(1 + 2 / sc) / gamma(1 + 1 / sc)^2 - 1)
    list(mean = dmean, sd = dmean * cv)
  } else {
    dmean <- exp(loc)
    list(mean = dmean, sd = dmean / sqrt(sc))
  }
}

#' Correct the ascertained CFR for outcome-dependent ascertainment
#'
#' Shifts the logit-scale CFR draws by `-log(r)`; `r` = 1 leaves the CFR
#' unchanged. See [summary.cfrnow_fit()] Details for the rationale.
#' @param logit_cfr Logit-scale CFR draws.
#' @param ascertainment_ratio Ratio `r` of the ascertainment probability of
#'   fatal to non-fatal cases.
#' @return CFR draws on the (0, 1) scale.
#' @noRd
.ascertainment_adjust <- function(logit_cfr, ascertainment_ratio) {
  stats::plogis(logit_cfr - log(ascertainment_ratio))
}

# Does the cfr sub-formula carry predictors beyond an intercept?
.cfr_is_grouped <- function(object) {
  cf <- object$formula$pforms$cfr
  if (is.null(cf)) {
    return(FALSE)
  }
  tt <- stats::terms(cf)
  length(attr(tt, "term.labels")) > 0L || attr(tt, "intercept") == 0L
}

# Per-group CFR draws for a `cfr ~ group` fit: one column per distinct
# combination of the (factor/character) cfr predictors, labelled by it. Predicts
# the cfr dpar on the fitted rows (parameterisation-agnostic, includes any
# group-level effects), then applies the ascertainment shift.
.cfr_group_draws <- function(object, ascertainment_ratio) {
  vars <- setdiff(all.vars(object$formula$pforms$cfr), "cfr")
  dat <- object$data
  discrete <- all(vapply(
    dat[vars],
    function(x) is.factor(x) || is.character(x) || is.logical(x),
    logical(1)
  ))
  if (!discrete) {
    stop("summary() reports per-group CFRs for factor or character `cfr` ",
      "predictors; for a continuous predictor use brms::posterior_epred() ",
      "with a newdata grid.",
      call. = FALSE
    )
  }
  ep <- brms::posterior_epred(object, dpar = "cfr")
  ep <- stats::plogis(stats::qlogis(ep) - log(ascertainment_ratio))
  key <- interaction(dat[vars], drop = TRUE, sep = ".")
  reps <- which(!duplicated(key))
  labs <- as.character(key)[reps]
  o <- order(labs)
  mat <- ep[, reps[o], drop = FALSE]
  colnames(mat) <- paste0("cfr[", labs[o], "]")
  mat
}

.cfr_quantities <- function(object, ascertainment_ratio = 1) {
  dr <- posterior::as_draws_df(object)
  grouped <- .cfr_is_grouped(object)
  if (!grouped && !"b_cfr_Intercept" %in% posterior::variables(dr)) {
    stop("summary() needs an intercept-only or `cfr ~ group` fit.",
      call. = FALSE
    )
  }
  fam <- object$cfrnow$family
  scale2 <- if (fam == "lognormal") "sigma" else "shape"
  d <- .delay_moments(
    dr, "b_Intercept",
    paste0("b_", scale2, "_Intercept"), fam
  )
  res <- if (grouped) {
    as.data.frame(
      .cfr_group_draws(object, ascertainment_ratio),
      check.names = FALSE
    )
  } else {
    data.frame(
      cfr = .ascertainment_adjust(dr[["b_cfr_Intercept"]], ascertainment_ratio)
    )
  }
  res$delay_mean <- d$mean
  res$delay_sd <- d$sd
  if (isTRUE(object$cfrnow$use_recovery)) {
    rfam <- object$cfrnow$recovery_family
    rscale2 <- if (rfam == "lognormal") "sigma" else "shape"
    r <- .delay_moments(
      dr, "b_rmu_Intercept",
      paste0("b_r", rscale2, "_Intercept"), rfam
    )
    res$recovery_mean <- r$mean
    res$recovery_sd <- r$sd
  }
  # A grouped fit's cfr columns come from posterior_epred, which returns draws
  # in as_draws_df order, so they line up with these chain/iteration ids.
  res$.chain <- dr$.chain
  res$.iteration <- dr$.iteration
  res$.draw <- dr$.draw
  posterior::as_draws_df(res)
}

#' Summarise a mixture-cure CFR fit
#'
#' Reports the corrected CFR and the onset-to-death delay (mean and sd, in days)
#' as posterior quantiles with convergence diagnostics (`rhat`, `ess_bulk`). The
#' naive `deaths / cases` ratio is returned as an attribute; in real time it
#' underestimates the corrected CFR because not every fatal case has died by the
#' cut-off.
#'
#' When few deaths have resolved (a young outbreak) the CFR is only weakly
#' identified and its posterior stays close to the prior. This is reported via
#' the `cfr_low_information` attribute: `TRUE` when the CFR posterior sd exceeds
#' `info_tol` times the prior sd.
#'
#' For a `cfr ~ group` fit the CFR varies by group, so one `cfr[<group>]` row is
#' reported per group (grouping predictors must be factors or characters), and
#' the `cfr_low_information` flag is `NA` (only defined for a single CFR).
#'
#' The CFR the model fits is the fatality risk among *ascertained* cases. When
#' ascertainment is outcome-dependent -- fatal and non-fatal cases entering the
#' line list at different rates -- this differs from the population CFR.
#' `ascertainment_ratio` (`r`) is the ratio of the ascertainment probability of
#' fatal to non-fatal cases; the reported CFR is shifted on the logit scale by
#' `-log(r)`, so `r` > 1 (fatal cases over-ascertained) lowers it and `r` < 1
#' (e.g. deaths not linked back to cases) raises it. It is supplied, not fitted,
#' and defaults to 1 (no correction); because the correction is a post-hoc logit
#' shift, sweep a range of `r` to show its leverage rather than trusting a
#' single value.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @param probs Quantiles to report.
#' @param info_tol Low-information threshold: flag when the CFR posterior sd is
#'   more than this fraction of the prior sd. Defaults to 0.9.
#' @param ascertainment_ratio Ratio `r` of the ascertainment probability of
#'   fatal to non-fatal cases (see Details). A single positive number; defaults
#'   to 1 (no correction).
#' @param ... Unused.
#' @return A data frame with one row per quantity: `cfr` (or one `cfr[<group>]`
#'   row per group for a `cfr ~ group` fit), `delay_mean` and `delay_sd`,
#'   carrying `naive_cfr`, `n_cases`, `n_deaths`, `cfr_prior_sd`,
#'   `cfr_low_information` and `ascertainment_ratio` attributes.
#' @family fit
#' @export
summary.cfrnow_fit <- function(object, probs = c(0.025, 0.5, 0.975),
                               info_tol = 0.9, ascertainment_ratio = 1, ...) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  if (!is.numeric(ascertainment_ratio) || length(ascertainment_ratio) != 1 ||
        !is.finite(ascertainment_ratio) || ascertainment_ratio <= 0) {
    stop("`ascertainment_ratio` must be a single positive number.",
      call. = FALSE
    )
  }
  qs <- .cfr_quantities(object, ascertainment_ratio)
  qcols <- paste0("q", probs * 100)
  sm <- posterior::summarise_draws(
    qs,
    mean = mean,
    stats::setNames(lapply(probs, function(p) {
      function(x) stats::quantile(x, p, names = FALSE)
    }), qcols),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )
  out <- as.data.frame(sm)
  names(out)[1] <- "quantity"

  # Weak identification is a single-CFR diagnostic: undo the shift and measure
  # the CFR spread on the fitted (r = 1) scale against the r = 1 prior sd.
  # It is not defined for a grouped fit, where the CFR varies by group.
  prior_sd <- object$cfrnow$cfr_prior_sd
  low_info <- NA
  if ("cfr" %in% posterior::variables(qs)) {
    cfr_obs <- stats::plogis(
      stats::qlogis(posterior::extract_variable(qs, "cfr")) +
        log(ascertainment_ratio)
    )
    low_info <- !is.na(prior_sd) && stats::sd(cfr_obs) > info_tol * prior_sd
  }
  attr(out, "naive_cfr") <- naive_cfr(
    object$cfrnow$n_deaths,
    object$cfrnow$n_cases
  )
  attr(out, "n_cases") <- object$cfrnow$n_cases
  attr(out, "n_deaths") <- object$cfrnow$n_deaths
  attr(out, "cfr_prior_sd") <- prior_sd
  attr(out, "cfr_low_information") <- low_info
  attr(out, "ascertainment_ratio") <- ascertainment_ratio
  out
}

#' @rdname summary.cfrnow_fit
#' @param x A `cfrnow_fit`.
#' @export
print.cfrnow_fit <- function(x, ...) {
  s <- summary(x)
  message("<cfrnow_fit> ", x$cfrnow$family, " delay")
  message(
    "  cases: ", x$cfrnow$n_cases,
    "   deaths by cut-off: ", x$cfrnow$n_deaths,
    "   naive CFR: ", round(attr(s, "naive_cfr"), 3)
  )
  print(s)
  if (isTRUE(attr(s, "cfr_low_information"))) {
    message(
      "  ! CFR only weakly identified (posterior close to prior); ",
      "few resolved deaths - interpret with caution."
    )
  }
  if (any(s$rhat > 1.01, na.rm = TRUE)) {
    message("  ! some Rhat > 1.01 - chains may not have converged.")
  }
  invisible(x)
}
