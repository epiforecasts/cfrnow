# Posterior-predictive check for a cfrnow fit

Draws replicate line-list outcomes from the posterior and compares them
with the data the model was fit to, replaying the real-time observation
process (the onset-to-event timing and the cut-off censoring). Two
checks are available: the count of observed deaths by the cut-off (plus
recoveries for a two-outcome fit), and the distribution of the observed
onset-to-death delays.

## Usage

``` r
pp_check_cfr(object, type = c("counts", "delay"), ndraws = 100)
```

## Arguments

- object:

  A `cfrnow_fit` from
  [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md).

- type:

  Which check to plot: `"counts"` (default) or `"delay"`.

- ndraws:

  Number of posterior draws to replicate over. Defaults to 100, or fewer
  if the fit has fewer draws.

## Value

A `ggplot` object.

## Details

The check reuses the fit's own posterior draws of the CFR and the delay,
so it works for covariate and time-varying `cfr ~ ...` fits as well as
intercept-only ones. It needs the observation cut-off, which
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md)
records when the data come from
[`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md);
a retrospective fit (`obs_time = NULL`) has no truncation to replay, so
every fatal case shows up as a death.

Only the quantities the model generates are checked: the observed death
(and recovery) counts and the observed onset-to-death delays. The split
of the remaining cases into censored versus untimed-resolved is not part
of the generative model, so it is left out.

The check is stochastic: it subsamples posterior draws and simulates
outcomes, so call [`set.seed()`](https://rdrr.io/r/base/Random.html)
beforehand for reproducible plots.

## See also

Other fit:
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md),
[`summary.cfrnow_fit()`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
fit <- fit_cfr(d, delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)))
pp_check_cfr(fit, type = "counts")
pp_check_cfr(fit, type = "delay")
} # }
```
