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

# Distribution function of a delay at times `q`, censored for the day of onset
# -- the censored branch of cure_lpmf_death.stan / cure_lpmf_two_outcome.stan,
# evaluated here with primarycensored's R censored-CDF rather than its Stan
# one. `loc` and `scale` are length-ndraws vectors of natural-scale delay
# parameters for a single case (a column of posterior_linpred(transform =
# TRUE)), and `pwindow` is that case's onset window. `delay_max` truncates the
# delay where the likelihood does (Inf when unbounded): the censored CDF is
# normalised by its value at the bound, as pprimarycensored() does for a
# finite `D`, which it will not do for vectors of parameters. `q` and
# `pwindow` have to match the parameter vectors' length for
# pprimarycensored() to vectorise correctly. Returns an ndraws x length(q)
# matrix.
.delay_cdf <- function(family, q, pwindow, loc, scale, delay_max = Inf) {
  nd <- length(loc)
  np <- .delay_native_params(family, loc, scale)
  cdf <- function(x) {
    do.call(primarycensored::pprimarycensored, c(
      list(
        q = rep(x, nd), pdist = np$pdist, pwindow = rep(pwindow, nd),
        L = -Inf, D = Inf, check = FALSE
      ),
      np$args
    ))
  }
  at_max <- if (is.finite(delay_max)) cdf(delay_max) else 1
  f <- vapply(q, function(x) {
    if (x >= delay_max) rep(1, nd) else cdf(x) / at_max
  }, numeric(nd))
  matrix(f, nd, length(q))
}

# Survivor function of a delay at each case's own follow-up time: `loc` and
# `scale` are ndraws x ncases matrices, and `y` and `pwindow` are
# length-ncases vectors, one value per case.
.delay_survivor <- function(family, y, pwindow, loc, scale, delay_max = Inf) {
  s <- vapply(seq_along(y), function(j) {
    1 - .delay_cdf(family, y[j], pwindow[j], loc[, j], scale[, j], delay_max)
  }, numeric(nrow(loc)))
  matrix(s, nrow(loc), length(y))
}

# Bayes update of `prob` on a case's follow-up: p * U_D(y) / (p * U_D(y) +
# (1 - p) * U_R(y)) for a still-unresolved (censored) case, else the
# deterministic outcome (1 for a death, 0 for a recovery or resolved
# non-death). U is the probability of still being unresolved on that outcome's
# path: the survivor function S, or l + (1 - l) S when a case can be lost with
# probability l, as in the fit's likelihood. `prob`, `s_d`, `s_r`, `loss_d` and
# `loss_r` are ndraws x ncases matrices (`s_r` may be the scalar 1 when there is
# no recovery delay, and either loss NULL when that outcome is never lost);
# `outcome` has length ncases.
.pi_death <- function(prob, s_d, s_r, outcome, loss_d = NULL, loss_r = NULL) {
  unresolved <- function(s, loss) {
    if (is.null(loss)) s else loss + (1 - loss) * s
  }
  u_d <- unresolved(s_d, loss_d)
  u_r <- unresolved(s_r, loss_r)
  pi_death <- prob * u_d / (prob * u_d + (1 - prob) * u_r)
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
#' delays, censored for the day of onset and truncated at each delay's `max`,
#' and `S_R = 1` when the fit has no `recovery_delay`. A fit with a
#' `loss_prior` replaces each survivor function `S` with `l + (1 - l) S`, where
#' `l` is the probability that a case with that outcome is lost to follow-up.
#' A case already resolved by the cut-off has a deterministic `pi`: 1 for an
#' observed death, 0 for a recovery or a resolved non-death.
#'
#' A case with little follow-up (small `y`) has `S_D(y)` close to 1, so its
#' `pi` sits close to its (covariate-specific) `prob`; only once it has been
#' followed for close to the typical delay does `pi` move away from `prob`
#' towards 0 or 1. Averaging `pi` over cases grouped by onset date therefore
#' gives the expected fatal fraction for each onset period without needing a
#' time trend in `prob`. For intervals on the fraction that actually dies, draw
#' each case's outcome from its `pi` before averaging. Estimates for the most
#' recent periods, whose cases have had little follow-up, lean towards the
#' fitted `prob` for those cases' covariates rather than their eventual
#' outcome.
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
#' # real-time CFR by onset week: draw each case's outcome from its pi, then
#' # take the fatal fraction of the cases in each week
#' pi <- posterior::as_draws_matrix(pi)
#' fatal <- matrix(rbinom(length(pi), 1, pi), nrow(pi))
#' week <- cut(fit$cfrnow$onset, "week")
#' by_week <- vapply(split(seq_along(week), week), function(i) {
#'   rowMeans(fatal[, i, drop = FALSE])
#' }, numeric(nrow(fatal)))
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
  s_d <- .delay_survivor(fam, d$y, d$pwindow, lp("mu"), lp(scale_dpar),
    delay_max = object$cfrnow$delay_max %||% Inf
  )

  s_r <- 1
  loss_r <- NULL
  if (isTRUE(object$cfrnow$use_recovery)) {
    rfam <- object$cfrnow$recovery_family
    rscale_dpar <- if (rfam == "lognormal") "rsigma" else "rshape"
    s_r <- .delay_survivor(rfam, d$y, d$pwindow, lp("rmu"), lp(rscale_dpar),
      delay_max = object$cfrnow$recovery_max %||% Inf
    )
    loss_r <- .ppc_loss(object, "recovery", lp, dim(prob))
  }

  pi_death <- .pi_death(prob, s_d, s_r, d$outcome,
    loss_d = .ppc_loss(object, "death", lp, dim(prob)),
    loss_r = loss_r
  )
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
