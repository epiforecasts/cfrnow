# Fit the real-time mixture-cure CFR model

Estimates a real-time case fatality ratio from line-list data. Each case
is fatal with probability `prob` and, when fatal, dies at an
onset-to-death `delay`; cases still alive at the cut-off are
right-censored, correcting the downward bias of the naive deaths / cases
ratio. The delay and the `prob` prior are given as distspec
distributions.

## Usage

``` r
fit_cfr(
  data,
  delay = LogNormal(meanlog = Normal(2, 1), sdlog = Normal(0.5, 0.3)),
  prob_prior = Beta(1, 1),
  recovery_delay = NULL,
  loss_prior = NULL,
  formula = mu ~ 1,
  ...,
  cfr_prior = lifecycle::deprecated()
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
  priors. A `max` (e.g.
  `LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15), max = 30)`) truncates
  the delay there, as it does everywhere else in distspec: recorded
  delays then run from 0 to `max - 1` days.

- prob_prior:

  Prior on `prob` as a
  [`distspec::Beta()`](https://epiforecasts.io/distspec/reference/Beta.html).
  Defaults to `Beta(1, 1)`.

- recovery_delay:

  Optional onset-to-recovery delay (same form as `delay`) for the
  two-outcome fit; may use a different family from `delay`.

- loss_prior:

  Optional prior on the probability that a case is lost to follow-up, so
  its outcome never reaches the line list. A
  [`distspec::Beta()`](https://epiforecasts.io/distspec/reference/Beta.html)
  gives one probability whatever the outcome; a list with `death` and
  `recovery` entries, each a
  [`Beta()`](https://epiforecasts.io/distspec/reference/Beta.html) or a
  fixed number, gives them their own (e.g.
  `list(death = 0, recovery = Beta(1, 1))` where every death is
  recorded). Needs recorded recoveries and a `recovery_delay`. `NULL`
  (the default) assumes every outcome is eventually recorded.

- formula:

  A `brms` formula for the delay location `mu` and, optionally, `prob`
  (`prob ~ ...`). Defaults to `mu ~ 1`. `prob_prior` normally lands on
  the `prob` intercept; when the `prob` formula drops the intercept
  (e.g. `prob ~ 0 + group`, one logit-probability per group) it is
  placed on those coefficients instead. It then applies to every `prob`
  coefficient, so an intercept-free formula should carry only factor
  terms. A `cfr ~ ...` sub-formula is still accepted and translated to
  `prob ~ ...` with a deprecation warning.

- ...:

  Passed to
  [`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html)
  and on to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  (e.g. `chains`, `iter`, `backend`, `seed`).

- cfr_prior:

  **\[deprecated\]** Use `prob_prior` instead.

## Value

A `brmsfit` with class `cfrnow_fit`; summarise with
[`summary()`](https://rdrr.io/r/base/summary.html).

## Details

The delay's native parameters may each be a fixed number (held fixed;
fixing the whole delay gives the Ghani/Nishiura estimator) or a
[`Normal()`](https://epiforecasts.io/distspec/reference/Normal.html)
prior (co-estimated). The family
([`LogNormal()`](https://epiforecasts.io/distspec/reference/LogNormal.html),
[`Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html) or
[`Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html))
sets the delay distribution. `prob_prior` is a
[`Beta()`](https://epiforecasts.io/distspec/reference/Beta.html); it
matters because `prob` is weakly identified early on (`Beta(1, 1)` is
uniform, `Beta(1, 9)` favours a low probability, `Beta(6.6, 13.4)` suits
a high-fatality pathogen).

`loss_prior` adds a probability that a case is lost to follow-up, so its
outcome never reaches the line list. A recorded death or recovery then
also says the case was kept, and a case still unresolved at the cut-off
was either lost or genuinely unresolved. Without it, a case unresolved
far longer than the delays allow has almost no probability, which pulls
the delay's tail out or stops the sampler from starting.

Whether loss can be estimated depends on how it relates to the outcome.
One probability for both (`loss_prior = Beta(1, 1)`) is identified by
the recorded deaths, the recorded recoveries and the cases that stay
unresolved, and it assumes a lost case's CFR matches everyone else's.
Fixing one half (`list(death = 0, recovery = Beta(1, 1))`, where a death
is always written down and a discharge may not be) is identified too.
Estimating both leaves one parameter unidentified: the priors decide the
split and `fit_cfr()` warns, so treat that fit as a sensitivity
analysis. A `loss_prior` needs recorded recoveries and a
`recovery_delay` to time them, because only timed recoveries separate
loss from the outcome probability; without them, `fit_cfr()` stops.

Where a case has a date it was last known unresolved, censoring it there
through
[`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)'s
`last_contact_date` uses that timing instead of inferring the loss.

A delay's `max` truncates it: the likelihood renormalises the
distribution over `[0, max)` days, matching how distspec discretises the
same object, so a delay recorded as `max` days or longer has no
probability. The estimated parameters still describe the untruncated
family, which is what [`summary()`](https://rdrr.io/r/base/summary.html)
reports as `delay_mean` and `delay_sd`. A recorded delay outside the
bound, or a case unresolved for as long as every bound that could still
apply to it, cannot come from such a model, so `fit_cfr()` stops rather
than letting the sampler fail.

The model is fitted through
[`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html),
so covariates (or a smooth time effect) can be put on `prob` or the
delay through `formula`, e.g. `formula = brms::bf(mu ~ 1, prob ~ age)`.
When the line list carries recovery dates, pass a `recovery_delay` to
fit the two-outcome model that also times recoveries.

`prob` is a case fatality ratio when the line list runs from onset to
death, as
[`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
expects; fed a line list with a different outcome (e.g.
hospitalisation), the same model fits that outcome's probability
instead.

## See also

Other fit:
[`posterior_prob_death()`](https://epiforecasts.io/cfrnow/reference/posterior_prob_death.md),
[`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md),
[`summary.cfrnow_fit()`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
otd <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))

# co-estimated delay
fit <- fit_cfr(d, delay = otd, prob_prior = Beta(1, 1), backend = "cmdstanr")
summary(fit)

# fixed delay (Ghani/Nishiura), a gamma family, and a prob covariate
fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51))
fit_cfr(d, delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)))
fit_cfr(d, delay = otd, formula = brms::bf(mu ~ 1, prob ~ group))
} # }
```
