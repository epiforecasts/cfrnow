# Fit the real-time mixture-cure CFR model

Estimates a real-time case fatality ratio from line-list data. Each case
is fatal with probability `cfr` and, when fatal, dies at an
onset-to-death `delay`; cases still alive at the cut-off are
right-censored, correcting the downward bias of the naive deaths / cases
ratio. The delay and the CFR prior are given as distspec distributions.

## Usage

``` r
fit_cfr(
  data,
  delay = LogNormal(meanlog = Normal(2, 1), sdlog = Normal(0.5, 0.3)),
  cfr_prior = Beta(1, 1),
  recovery_delay = NULL,
  formula = mu ~ 1,
  ...
)
```

## Arguments

- data:

  Output of
  [`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md),
  or an `epidist_cure_model` / data frame with `y`, `outcome`,
  `pwindow`, `swindow`.

- delay:

  Onset-to-death delay as a distspec distribution
  ([`distspec::LogNormal()`](https://epiforecasts.io/distspec/reference/LogNormal.html),
  [`distspec::Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html)
  or
  [`distspec::Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html))
  whose native parameters are fixed numbers or
  [`Normal()`](https://epiforecasts.io/distspec/reference/Normal.html)
  priors.

- cfr_prior:

  CFR prior as a
  [`distspec::Beta()`](https://epiforecasts.io/distspec/reference/Beta.html).
  Defaults to `Beta(1, 1)`.

- recovery_delay:

  Optional onset-to-recovery delay (same form as `delay`) for the
  two-outcome fit; may use a different family from `delay`.

- formula:

  A `brms` formula for the delay location `mu` and, optionally, the CFR
  (`cfr ~ ...`). Defaults to `mu ~ 1`. The `cfr_prior` normally lands on
  the CFR intercept; when the CFR formula drops the intercept (e.g.
  `cfr ~ 0 + group`, one logit-CFR per group) it is placed on those
  coefficients instead. It then applies to every `cfr` coefficient, so
  an intercept-free formula should carry only factor terms.

- ...:

  Passed to
  [`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html)
  and on to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  (e.g. `chains`, `iter`, `backend`, `seed`).

## Value

A `brmsfit` with class `cfrnow_fit`; summarise with
[`cmdstanr::summary()`](https://mc-stan.org/cmdstanr/reference/fit-method-summary.html).

## Details

The delay's native parameters may each be a fixed number (held fixed;
fixing the whole delay gives the Ghani/Nishiura estimator) or a
[`Normal()`](https://epiforecasts.io/distspec/reference/Normal.html)
prior (co-estimated). The family
([`LogNormal()`](https://epiforecasts.io/distspec/reference/LogNormal.html),
[`Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html) or
[`Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html))
sets the delay distribution. `cfr_prior` is a
[`Beta()`](https://epiforecasts.io/distspec/reference/Beta.html); it
matters because the CFR is weakly identified early on (`Beta(1, 1)` is
uniform, `Beta(1, 9)` favours a low CFR, `Beta(6.6, 13.4)` suits a
high-fatality pathogen).

The model is fitted through
[`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html),
so covariates (or a smooth time effect) can be put on the CFR or the
delay through `formula`, e.g. `formula = brms::bf(mu ~ 1, cfr ~ age)`.
When the line list carries recovery dates, pass a `recovery_delay` to
fit the two-outcome model that also times recoveries.

## See also

Other fit:
[`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md),
[`summary.cfrnow_fit()`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
otd <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))

# co-estimated delay
fit <- fit_cfr(d, delay = otd, cfr_prior = Beta(1, 1), backend = "cmdstanr")
summary(fit)

# fixed delay (Ghani/Nishiura), a gamma family, and a CFR covariate
fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51))
fit_cfr(d, delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)))
fit_cfr(d, delay = otd, formula = brms::bf(mu ~ 1, cfr ~ group))
} # }
```
