# Native distribution parameters for a delay family, from posterior_linpred's
# natural-scale `mu` and second parameter (the same reparameterisation
# .family_stan_param() reads out of brms for cure_model.R's Stan template):
# lognormal is already (meanlog, sdlog); gamma and weibull are mean-
# parameterised in brms, so are converted to (shape, rate) and (shape, scale).
.delay_native_params <- function(family, loc, scale) {
  if (family == "lognormal") {
    list(pdist = stats::plnorm, args = list(meanlog = loc, sdlog = scale))
  } else if (family == "weibull") {
    list(
      pdist = stats::pweibull,
      args = list(shape = scale, scale = loc / gamma(1 + 1 / scale))
    )
  } else {
    list(pdist = stats::pgamma, args = list(shape = scale, rate = scale / loc))
  }
}

# Survivor function of a delay, censored for the day of onset -- the censored
# branch of cure_lpmf_death.stan / cure_lpmf_two_outcome.stan, evaluated here
# with primarycensored's R censored-CDF rather than its Stan one. `loc` and
# `scale` are ndraws x ncases matrices of natural-scale delay parameters (as
# returned by posterior_linpred(transform = TRUE)); `y` and `pwindow` are
# length-ncases vectors, one value per case. `q` and `pwindow` have to match
# the parameter vectors' length for pprimarycensored() to vectorise correctly.
.delay_survivor <- function(family, y, pwindow, loc, scale) {
  nd <- nrow(loc)
  s <- matrix(NA_real_, nd, ncol(loc))
  for (j in seq_len(ncol(loc))) {
    np <- .delay_native_params(family, loc[, j], scale[, j])
    fbar <- do.call(primarycensored::pprimarycensored, c(
      list(
        q = rep(y[j], nd), pdist = np$pdist, pwindow = rep(pwindow[j], nd),
        L = -Inf, D = Inf, check = FALSE
      ),
      np$args
    ))
    s[, j] <- 1 - fbar
  }
  s
}

# Bayes update of `prob` on a case's follow-up: p * S_D(y) / (p * S_D(y) +
# (1 - p) * S_R(y)) for a still-unresolved (censored) case, else the
# deterministic outcome (1 for a death, 0 for a recovery or resolved
# non-death). `prob`, `s_d` and `s_r` are ndraws x ncases matrices (`s_r` may
# be the scalar 1 when there is no recovery delay); `outcome` has length ncases.
.pi_death <- function(prob, s_d, s_r, outcome) {
  pi_death <- prob * s_d / (prob * s_d + (1 - prob) * s_r)
  pi_death[, outcome == .CURE_DEATH] <- 1
  pi_death[, outcome == .CURE_RECOVERY | outcome == .CURE_RESOLVED] <- 0
  pi_death
}

#' Posterior probability of death for each case
#'
#' [brms::posterior_linpred()] gives each case's `prob`, the probability of
#' death before anything about the case's outcome has been observed. For a
#' case still unresolved at the cut-off after `y` days of follow-up, this
#' updates `prob` with that follow-up by Bayes' rule (the censored branch of
#' the cure likelihood):
#' \deqn{\pi = \frac{p \, S_D(y)}{p \, S_D(y) + (1 - p) \, S_R(y)}}
#' where `S_D` and `S_R` are the survivor functions of the death and recovery
#' delays, censored for the day of onset, and `S_R = 1` when the fit has no
#' `recovery_delay`. A case already resolved by the cut-off has a
#' deterministic `pi`: 1 for an observed death, 0 for a recovery or a resolved
#' non-death.
#'
#' A case with little follow-up (small `y`) has `S_D(y)` close to 1, so its
#' `pi` sits close to its (covariate-specific) `prob`; only once it has been
#' followed for close to the typical delay does `pi` move away from `prob`
#' towards 0 or 1. Averaging `pi` over cases grouped by onset date therefore
#' gives a real-time CFR for each onset period, with credible intervals and
#' without needing a time trend in `prob` -- but estimates for the most recent
#' periods, whose cases have had little follow-up, lean towards the fitted
#' `prob` for those cases' covariates rather than their eventual outcome.
#'
#' Works for any of the supported delay families and any `prob` or delay
#' `formula`, since the per-case parameters come from `posterior_linpred()`.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @return A `draws_df` with one `pi[<i>]` column per case, in the row order
#'   of `object$data`, holding posterior draws of the probability of death.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' fit <- fit_cfr(d, delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)))
#' pi <- posterior_prob_death(fit)
#'
#' # real-time CFR by onset week: average pi over the cases in each week
#' week <- cut(fit$cfrnow$onset, "week")
#' by_week <- vapply(split(seq_along(week), week), function(i) {
#'   rowMeans(pi[, i, drop = FALSE])
#' }, numeric(nrow(pi)))
#' apply(by_week, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
#' }
#' @family fit
#' @export
posterior_prob_death <- function(object) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  d <- object$data
  fam <- object$cfrnow$family
  scale_dpar <- if (fam == "lognormal") "sigma" else "shape"
  lp <- function(dpar) {
    brms::posterior_linpred(object, dpar = dpar, transform = TRUE)
  }

  prob <- lp("prob")
  s_d <- .delay_survivor(fam, d$y, d$pwindow, lp("mu"), lp(scale_dpar))

  s_r <- 1
  if (isTRUE(object$cfrnow$use_recovery)) {
    rfam <- object$cfrnow$recovery_family
    rscale_dpar <- if (rfam == "lognormal") "rsigma" else "rshape"
    s_r <- .delay_survivor(rfam, d$y, d$pwindow, lp("rmu"), lp(rscale_dpar))
  }

  pi_death <- .pi_death(prob, s_d, s_r, d$outcome)
  colnames(pi_death) <- paste0("pi[", seq_len(ncol(pi_death)), "]")

  # Draws come out in as_draws_df() row order (see .cfr_quantities()), so the
  # chain/iteration/draw ids line up with the pi columns computed above.
  dr <- posterior::as_draws_df(object)
  res <- as.data.frame(pi_death, check.names = FALSE)
  res$.chain <- dr$.chain
  res$.iteration <- dr$.iteration
  res$.draw <- dr$.draw
  posterior::as_draws_df(res)
}
