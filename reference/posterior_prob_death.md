# Posterior probability of death for each case

[`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)
gives each case's `prob`, the probability of death before anything about
the case's outcome has been observed. For a case still unresolved at the
cut-off after `y` days of follow-up, this updates `prob` with that
follow-up by Bayes' rule (the censored branch of the cure likelihood):
\$\$\pi = \frac{p \\ S_D(y)}{p \\ S_D(y) + (1 - p) \\ S_R(y)}\$\$ where
`S_D` and `S_R` are the survivor functions of the death and recovery
delays, censored for the day of onset and truncated at each delay's
`max`, and `S_R = 1` when the fit has no `recovery_delay`. A fit with a
`loss_prior` replaces each survivor function `S` with `l + (1 - l) S`,
where `l` is the probability that a case with that outcome is lost to
follow-up. A case already resolved by the cut-off has a deterministic
`pi`: 1 for an observed death, 0 for a recovery or a resolved non-death.

## Usage

``` r
posterior_prob_death(object)
```

## Arguments

- object:

  A `cfrnow_fit` from
  [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md).

## Value

A `draws_df` with one `pi[<i>]` column per case, in the row order of
`object$data`, holding posterior draws of the probability of death.

## Details

A case with little follow-up (small `y`) has `S_D(y)` close to 1, so its
`pi` sits close to its (covariate-specific) `prob`; only once it has
been followed for close to the typical delay does `pi` move away from
`prob` towards 0 or 1. Averaging `pi` over cases grouped by onset date
therefore gives the expected fatal fraction for each onset period
without needing a time trend in `prob`. For intervals on the fraction
that actually dies, draw each case's outcome from its `pi` before
averaging. Estimates for the most recent periods, whose cases have had
little follow-up, lean towards the fitted `prob` for those cases'
covariates rather than their eventual outcome.

Works for any of the supported delay families and any `prob` or delay
`formula`, since the per-case parameters come from
`posterior_linpred()`.

## See also

Other fit:
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md),
[`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md),
[`summary.cfrnow_fit()`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
fit <- fit_cfr(d, delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)))
pi <- posterior_prob_death(fit)

# real-time CFR by onset week: draw each case's outcome from its pi, then
# take the fatal fraction of the cases in each week
pi <- posterior::as_draws_matrix(pi)
fatal <- matrix(rbinom(length(pi), 1, pi), nrow(pi))
week <- cut(fit$cfrnow$onset, "week")
by_week <- vapply(split(seq_along(week), week), function(i) {
  rowMeans(fatal[, i, drop = FALSE])
}, numeric(nrow(fatal)))
apply(by_week, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
} # }
```
