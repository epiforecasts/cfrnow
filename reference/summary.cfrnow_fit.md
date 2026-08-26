# Summarise a mixture-cure CFR fit

Reports the corrected CFR and the onset-to-death delay (mean and sd, in
days) as posterior quantiles with convergence diagnostics (`rhat`,
`ess_bulk`). The naive `deaths / cases` ratio is returned as an
attribute; in real time it underestimates the corrected CFR because not
every fatal case has died by the cut-off.

## Usage

``` r
# S3 method for class 'cfrnow_fit'
summary(
  object,
  probs = c(0.025, 0.5, 0.975),
  info_tol = 0.9,
  ascertainment_ratio = 1,
  ...
)

# S3 method for class 'cfrnow_fit'
print(x, ...)
```

## Arguments

- object:

  A `cfrnow_fit` from
  [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md).

- probs:

  Quantiles to report.

- info_tol:

  Low-information threshold: flag when the CFR posterior sd is more than
  this fraction of the prior sd. Defaults to 0.9.

- ascertainment_ratio:

  Ratio `r` of the ascertainment probability of fatal to non-fatal cases
  (see Details). A single positive number; defaults to 1 (no
  correction).

- ...:

  Unused.

- x:

  A `cfrnow_fit`.

## Value

A data frame with one row per quantity: `cfr` (or one `cfr[<group>]` row
per group for a `cfr ~ group` fit), `delay_mean` and `delay_sd`,
carrying `naive_cfr`, `n_cases`, `n_deaths`, `cfr_prior_sd`,
`cfr_low_information` and `ascertainment_ratio` attributes.

## Details

When few deaths have resolved (a young outbreak) the CFR is only weakly
identified and its posterior stays close to the prior. This is reported
via the `cfr_low_information` attribute: `TRUE` when the CFR posterior
sd exceeds `info_tol` times the prior sd.

For a `cfr ~ group` fit the CFR varies by group, so one `cfr[<group>]`
row is reported per group (grouping predictors must be factors or
characters), and the `cfr_low_information` flag is `NA` (only defined
for a single CFR).

The CFR the model fits is the fatality risk among *ascertained* cases.
When ascertainment is outcome-dependent – fatal and non-fatal cases
entering the line list at different rates – this differs from the
population CFR. `ascertainment_ratio` (`r`) is the ratio of the
ascertainment probability of fatal to non-fatal cases; the reported CFR
is shifted on the logit scale by `-log(r)`, so `r` \> 1 (fatal cases
over-ascertained) lowers it and `r` \< 1 (e.g. deaths not linked back to
cases) raises it. It is supplied, not fitted, and defaults to 1 (no
correction); because the correction is a post-hoc logit shift, sweep a
range of `r` to show its leverage rather than trusting a single value.

## See also

Other fit:
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md),
[`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md)
