# Summarise a mixture-cure CFR fit

Pools the posterior draws of the corrected CFR and the onset-to-death
delay (mean and standard deviation, in days) and reports quantile
summaries with convergence diagnostics (`rhat`, `ess_bulk`). The naive
`deaths / cases` ratio is returned as an attribute for comparison; in
real time it underestimates the corrected CFR because not every fatal
case has died by the cut-off. The `delay_mean`/`delay_sd` summaries
(days) are recovered from the native parameters as generated quantities,
so they are the same whatever family the delay used; with a fixed delay
they are constant.

## Usage

``` r
summarise_cfr(object, probs = c(0.025, 0.5, 0.975), info_tol = 0.9, ...)
```

## Arguments

- object:

  A `cfrnow_fit` from
  [`fit_cfr()`](https://sbfnk.github.io/cfrnow/reference/fit_cfr.md).

- probs:

  Quantiles to report.

- info_tol:

  Low-information threshold: flag when the CFR posterior sd is more than
  this fraction of the prior sd. Defaults to 0.9.

- ...:

  Unused.

## Value

A data frame with one row per summarised quantity (`cfr`, `delay_mean`,
`delay_sd`), carrying `naive_cfr`, `n_cases`, `n_deaths`, `cfr_prior_sd`
and `cfr_low_information` attributes.

## Details

When few deaths have resolved (a young outbreak), the CFR is only weakly
identified and its posterior stays close to the prior. This is reported
via the `cfr_low_information` attribute: `TRUE` when the CFR posterior
sd exceeds `info_tol` times the prior sd.

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(delay = dist.spec::LogNormal(2.4, 0.5))
fit <- fit_cfr(prepare_cfr_data(ll),
               delay = dist.spec::LogNormal(2.4, 0.5))
summarise_cfr(fit)
} # }
```
