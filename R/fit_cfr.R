#' Fit the real-time mixture-cure CFR model
#'
#' Estimates a real-time case fatality ratio from line-list data. Each case is
#' fatal with probability `prob` and, when fatal, dies at an onset-to-death
#' `delay`; cases still alive at the cut-off are right-censored, correcting the
#' downward bias of the naive deaths / cases ratio. The delay and the `prob`
#' prior are given as \pkg{distspec} distributions.
#'
#' The delay's native parameters may each be a fixed number (held fixed; fixing
#' the whole delay gives the Ghani/Nishiura estimator) or a `Normal()` prior
#' (co-estimated). The family (`LogNormal()`, `Gamma()` or `Weibull()`) sets the
#' delay distribution. `prob_prior` is a `Beta()`; it matters because `prob` is
#' weakly identified early on (`Beta(1, 1)` is uniform, `Beta(1, 9)` favours a
#' low probability, `Beta(6.6, 13.4)` suits a high-fatality pathogen).
#'
#' A delay's `max` truncates it: the likelihood renormalises the distribution
#' over `[0, max)` days, matching how \pkg{distspec} discretises the same
#' object, so a delay recorded as `max` days or longer has no probability. The
#' estimated parameters still describe the untruncated family, which is what
#' `summary()` reports as `delay_mean` and `delay_sd`. A recorded delay outside
#' the bound, or a case unresolved for as long as every bound that could still
#' apply to it, cannot come from such a model, so `fit_cfr()` stops rather than
#' letting the sampler fail.
#'
#' The model is fitted through [epidist::epidist()], so covariates (or a smooth
#' time effect) can be put on `prob` or the delay through `formula`, e.g.
#' `formula = brms::bf(mu ~ 1, prob ~ age)`. When the line list carries recovery
#' dates, pass a `recovery_delay` to fit the two-outcome model that also times
#' recoveries.
#'
#' `prob` is a case fatality ratio when the line list runs from onset to death,
#' as [prepare_cfr_data()] expects; fed a line list with a different outcome
#' (e.g. hospitalisation), the same model fits that outcome's probability
#' instead.
#'
#' @param data Output of [prepare_cfr_data()], or an `epidist_cure_model` /
#'   data frame with `y`, `outcome`, `pwindow`, `swindow`.
#' @param delay Onset-to-death delay as a \pkg{distspec} distribution
#'   ([distspec::LogNormal()], [distspec::Gamma()] or [distspec::Weibull()])
#'   whose native parameters are fixed numbers or `Normal()` priors. A `max`
#'   (e.g. `LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15), max = 30)`) truncates
#'   the delay there, as it does everywhere else in \pkg{distspec}: recorded
#'   delays then run from 0 to `max - 1` days.
#' @param prob_prior Prior on `prob` as a [distspec::Beta()]. Defaults to
#'   `Beta(1, 1)`.
#' @param recovery_delay Optional onset-to-recovery delay (same form as `delay`)
#'   for the two-outcome fit; may use a different family from `delay`.
#' @param formula A `brms` formula for the delay location `mu` and, optionally,
#'   `prob` (`prob ~ ...`). Defaults to `mu ~ 1`. `prob_prior` normally lands on
#'   the `prob` intercept; when the `prob` formula drops the intercept (e.g.
#'   `prob ~ 0 + group`, one logit-probability per group) it is placed on those
#'   coefficients instead. It then applies to every `prob` coefficient, so an
#'   intercept-free formula should carry only factor terms. A `cfr ~ ...`
#'   sub-formula is still accepted and translated to `prob ~ ...` with a
#'   deprecation warning.
#' @param ... Passed to [epidist::epidist()] and on to [brms::brm()]
#'   (e.g. `chains`, `iter`, `backend`, `seed`).
#' @param cfr_prior `r lifecycle::badge("deprecated")` Use `prob_prior` instead.
#' @return A `brmsfit` with class `cfrnow_fit`; summarise with [summary()].
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' otd <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))
#'
#' # co-estimated delay
#' fit <- fit_cfr(d, delay = otd, prob_prior = Beta(1, 1), backend = "cmdstanr")
#' summary(fit)
#'
#' # fixed delay (Ghani/Nishiura), a gamma family, and a prob covariate
#' fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51))
#' fit_cfr(d, delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)))
#' fit_cfr(d, delay = otd, formula = brms::bf(mu ~ 1, prob ~ group))
#' }
#' @family fit
#' @export
fit_cfr <- function(data,
                    delay = LogNormal(
                      meanlog = Normal(2, 1),
                      sdlog = Normal(0.5, 0.3)
                    ),
                    prob_prior = Beta(1, 1), recovery_delay = NULL,
                    formula = mu ~ 1, ...,
                    cfr_prior = lifecycle::deprecated()) {
  if (lifecycle::is_present(cfr_prior)) {
    lifecycle::deprecate_soft(
      "0.3.0", "fit_cfr(cfr_prior = )", "fit_cfr(prob_prior = )"
    )
    prob_prior <- cfr_prior
  }
  if (.formula_has_cfr(formula)) {
    lifecycle::deprecate_soft(
      "0.3.0", "fit_cfr(formula = 'must use `prob ~ ...`, not `cfr ~ ...`')"
    )
    formula <- .rename_cfr_formula(formula)
  }
  # Kept so posterior-predictive checks can replay the real-time truncation;
  # NA (retrospective) means every case is fully followed up.
  obs_time <- if (inherits(data, "cfrnow_data")) data$obs_time else as.Date(NA)
  cure <- as_epidist_cure_model(data)
  dd <- .delay_family_prior(delay, main = TRUE)
  dfam <- dd$family
  attr(cure, "delay_max") <- dd$max
  prob_class <- if (.prob_has_intercept(formula)) "Intercept" else "b"
  prior <- c(.prob_prior_to_brms(prob_prior, prob_class), dd$prior)
  rfam <- dfam
  if (isTRUE(attr(cure, "use_recovery"))) {
    if (is.null(recovery_delay)) {
      # no recovery delay supplied: treat recoveries as untimed resolutions
      cure$outcome[cure$outcome == .CURE_RECOVERY] <- .CURE_RESOLVED
      attr(cure, "use_recovery") <- FALSE
    } else {
      rd <- .delay_family_prior(recovery_delay, main = FALSE)
      rfam <- rd$family
      prior <- c(prior, rd$prior)
      attr(cure, "recovery_family") <- brms:::validate_family(rfam) # nolint
      attr(cure, "recovery_max") <- rd$max
    }
  }
  use_recovery <- isTRUE(attr(cure, "use_recovery"))
  .assert_within_max(cure, dd$max, if (use_recovery) attr(cure, "recovery_max"))
  fit <- epidist::epidist(cure,
    formula = formula, family = dfam,
    prior = prior, merge_priors = FALSE, ...
  )
  # brms drops rows with a missing formula covariate before fitting, so the
  # stored onset dates (used by pp_check_cfr() to replay real-time truncation)
  # must be subset to the rows it actually kept; brms preserves the row names
  # of `cure`, including after input subsetting. Warn too, since dropping cases
  # changes the estimand.
  used_rows <- match(rownames(fit$data), rownames(cure))
  n_dropped <- nrow(cure) - length(used_rows)
  if (n_dropped > 0) {
    warning(n_dropped, " case(s) dropped because a formula covariate is ",
      "missing (NA); the model was fitted on the remaining complete cases.",
      call. = FALSE
    )
  }
  fit$cfrnow <- list(
    n_cases = nrow(cure),
    n_deaths = sum(cure$outcome == .CURE_DEATH),
    family = brms:::validate_family(dfam)$family, # nolint
    use_recovery = use_recovery,
    recovery_family = if (use_recovery) {
      brms:::validate_family(rfam)$family # nolint
    } else {
      NA_character_
    },
    prob_prior_sd = .prob_prior_sd(prior),
    delay_max = dd$max,
    recovery_max = if (use_recovery) attr(cure, "recovery_max") else Inf,
    obs_time = obs_time,
    onset = if ("onset" %in% names(cure)) cure$onset[used_rows] else NULL
  )
  class(fit) <- c("cfrnow_fit", class(fit))
  fit
}

#' Prior sd of `prob` implied by a logit-scale Normal prior on its intercept
#'
#' Simulates the (0, 1)-scale sd from the `prob` intercept prior so [summary()]
#' can flag when the posterior is barely tighter than the prior. Returns `NA`
#' when the prior is not a recognised `normal(mean, sd)` on the prob intercept.
#' @param prior A `brmsprior`.
#' @return A numeric prior sd on the (0, 1) scale, or `NA`.
#' @noRd
.prob_prior_sd <- function(prior) {
  rows <- which(prior$dpar == "prob" & prior$class == "Intercept")
  if (length(rows) != 1) {
    return(NA_real_)
  }
  m <- regmatches(
    prior$prior[rows],
    regexec("normal\\(([^,]+),\\s*([^)]+)\\)", prior$prior[rows])
  )[[1]]
  if (length(m) != 3) {
    return(NA_real_)
  }
  draws <- stats::rnorm(1e5, as.numeric(m[2]), as.numeric(m[3]))
  stats::sd(stats::plogis(draws))
}

#' Check the data against the delays' upper bounds
#'
#' A bounded delay gives zero probability to anything longer than its `max`, so
#' a recorded delay past the bound, or a case still unresolved past every bound
#' that applies to it, cannot have come from the model and would make the fit
#' fail inside Stan. Stop with a message naming the cases instead.
#' @param cure An `epidist_cure_model`.
#' @param delay_max,recovery_max Upper bounds of the death and recovery delays;
#'   `Inf` (or `NULL`) when unbounded.
#' @return `TRUE`, invisibly.
#' @noRd
.assert_within_max <- function(cure, delay_max, recovery_max = NULL) {
  recovery_max <- recovery_max %||% Inf
  # primarycensored rejects a delay whose secondary window closes past the
  # bound, so the usable range is y + swindow <= max
  too_long <- function(code, mx) {
    sum(cure$outcome == code & cure$y + cure$swindow > mx)
  }
  n_death <- too_long(.CURE_DEATH, delay_max)
  if (n_death > 0) {
    stop(n_death, " death(s) with an onset-to-death delay that the delay's ",
      "max (", delay_max, " days) gives no probability. Raise the max, or ",
      "drop those records as data errors.",
      call. = FALSE
    )
  }
  n_recovery <- too_long(.CURE_RECOVERY, recovery_max)
  if (n_recovery > 0) {
    stop(n_recovery, " recovery(ies) with an onset-to-recovery delay that ",
      "the recovery delay's max (", recovery_max, " days) gives no ",
      "probability. Raise the max, or drop those records as data errors.",
      call. = FALSE
    )
  }
  # A censored case must still be able to resolve: it needs at least one of the
  # outcomes it could still have to remain possible beyond its follow-up.
  unresolved_max <- if (isTRUE(attr(cure, "use_recovery"))) {
    max(delay_max, recovery_max)
  } else {
    # death-only: an unresolved case may always be a survivor
    Inf
  }
  n_cens <- sum(cure$outcome == .CURE_CENSORED & cure$y >= unresolved_max)
  if (n_cens > 0) {
    stop(n_cens, " case(s) still unresolved at or past the delays' max (",
      unresolved_max, " days), which the model gives zero probability. ",
      "Drop them, or record them as resolved non-deaths.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
