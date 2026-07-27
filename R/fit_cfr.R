#' Fit the real-time mixture-cure CFR model
#'
#' Estimates a real-time case fatality ratio from line-list data. Each case is
#' fatal with probability `cfr` and, when fatal, dies at an onset-to-death
#' `delay`; cases still alive at the cut-off are right-censored, correcting the
#' downward bias of the naive deaths / cases ratio. The delay and the CFR prior
#' are given as \pkg{distspec} distributions.
#'
#' The delay's native parameters may each be a fixed number (held fixed; fixing
#' the whole delay gives the Ghani/Nishiura estimator) or a `Normal()` prior
#' (co-estimated). The family (`LogNormal()`, `Gamma()` or `Weibull()`) sets the
#' delay distribution. `cfr_prior` is a `Beta()`; it matters because the CFR is
#' weakly
#' identified early on (`Beta(1, 1)` is uniform, `Beta(1, 9)` favours a low CFR,
#' `Beta(6.6, 13.4)` suits a high-fatality pathogen).
#'
#' The model is fitted through [epidist::epidist()], so covariates (or a smooth
#' time effect) can be put on the CFR or the delay through `formula`, e.g.
#' `formula = brms::bf(mu ~ 1, cfr ~ age)`. When the line list carries recovery
#' dates, pass a `recovery_delay` to fit the two-outcome model that also times
#' recoveries.
#'
#' @param data Output of [prepare_cfr_data()], or an `epidist_cure_model` /
#'   data frame with `y`, `outcome`, `pwindow`, `swindow`.
#' @param delay Onset-to-death delay as a \pkg{distspec} distribution
#'   ([distspec::LogNormal()], [distspec::Gamma()] or [distspec::Weibull()])
#'   whose native parameters are fixed numbers or `Normal()` priors.
#' @param cfr_prior CFR prior as a [distspec::Beta()]. Defaults to `Beta(1, 1)`.
#' @param recovery_delay Optional onset-to-recovery delay (same form as `delay`)
#'   for the two-outcome fit; may use a different family from `delay`.
#' @param formula A `brms` formula for the delay location `mu` and, optionally,
#'   the CFR (`cfr ~ ...`). Defaults to `mu ~ 1`. The `cfr_prior` normally lands
#'   on the CFR intercept; when the CFR formula drops the intercept (e.g.
#'   `cfr ~ 0 + group`, one logit-CFR per group) it is placed on those
#'   coefficients instead. It then applies to every `cfr` coefficient, so an
#'   intercept-free formula should carry only factor terms.
#' @param ... Passed to [epidist::epidist()] and on to [brms::brm()]
#'   (e.g. `chains`, `iter`, `backend`, `seed`).
#' @return A `brmsfit` with class `cfrnow_fit`; summarise with [summary()].
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' otd <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))
#'
#' # co-estimated delay
#' fit <- fit_cfr(d, delay = otd, cfr_prior = Beta(1, 1), backend = "cmdstanr")
#' summary(fit)
#'
#' # fixed delay (Ghani/Nishiura), a gamma family, and a CFR covariate
#' fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51))
#' fit_cfr(d, delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)))
#' fit_cfr(d, delay = otd, formula = brms::bf(mu ~ 1, cfr ~ group))
#' }
#' @family fit
#' @export
fit_cfr <- function(data,
                    delay = LogNormal(
                      meanlog = Normal(2, 1),
                      sdlog = Normal(0.5, 0.3)
                    ),
                    cfr_prior = Beta(1, 1), recovery_delay = NULL,
                    formula = mu ~ 1, ...) {
  # Kept so posterior-predictive checks can replay the real-time truncation;
  # NA (retrospective) means every case is fully followed up.
  obs_time <- if (inherits(data, "cfrnow_data")) data$obs_time else as.Date(NA)
  cure <- as_epidist_cure_model(data)
  dd <- .delay_family_prior(delay, main = TRUE)
  dfam <- dd$family
  cfr_class <- if (.cfr_has_intercept(formula)) "Intercept" else "b"
  prior <- c(.cfr_prior_to_brms(cfr_prior, cfr_class), dd$prior)
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
    }
  }
  use_recovery <- isTRUE(attr(cure, "use_recovery"))
  fit <- epidist::epidist(cure,
    formula = formula, family = dfam,
    prior = prior, merge_priors = FALSE, ...
  )
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
    cfr_prior_sd = .cfr_prior_sd(prior),
    obs_time = obs_time,
    onset = if ("onset" %in% names(cure)) cure$onset else NULL
  )
  class(fit) <- c("cfrnow_fit", class(fit))
  fit
}

#' Prior sd of the CFR implied by a logit-scale Normal prior on its intercept
#'
#' Simulates the (0, 1)-scale sd from the `cfr` intercept prior so [summary()]
#' can flag when the posterior is barely tighter than the prior. Returns `NA`
#' when the prior is not a recognised `normal(mean, sd)` on the cfr intercept.
#' @param prior A `brmsprior`.
#' @return A numeric prior sd on the CFR scale, or `NA`.
#' @noRd
.cfr_prior_sd <- function(prior) {
  rows <- which(prior$dpar == "cfr" & prior$class == "Intercept")
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
