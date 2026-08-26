# Simulate a line list for testing and examples

Draws onset dates over a window, marks each case fatal with probability
`cfr`, and gives fatal cases an onset-to-death delay drawn from `delay`.
Pass a `recovery` delay to also simulate onset-to-recovery times for the
non-fatal cases and add a `recovery_date` column. Delays are distspec
distributions with fixed parameters, matching what
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md)
takes; you can specify them by mean and sd (e.g.
`LogNormal(mean = 12.75, sd = 7)`). Returns a line list with
`onset_date`, `death_date` (`NA` for non-fatal cases) and, when
`recovery` is given, `recovery_date` (`NA` for fatal cases). The full,
untruncated outcomes are simulated; pass the result to
[`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
with an `obs_time` to induce the real-time truncation.

## Usage

``` r
simulate_linelist(
  n = 200,
  cfr = 0.5,
  delay,
  recovery = NULL,
  onset_start = as.Date("2026-01-01"),
  onset_days = 60
)
```

## Arguments

- n:

  Number of cases.

- cfr:

  True case fatality ratio.

- delay:

  Onset-to-death delay: a distspec distribution
  ([`distspec::LogNormal()`](https://epiforecasts.io/distspec/reference/LogNormal.html),
  [`distspec::Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html)
  or
  [`distspec::Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html))
  with fixed parameters.

- recovery:

  Optional onset-to-recovery delay (same form as `delay`); when given,
  non-fatal cases get a `recovery_date`.

- onset_start:

  First possible onset date.

- onset_days:

  Width of the onset window (days); onsets are uniform over it.

## Value

A data frame with `onset_date`, `death_date` and, if `recovery` is
given, `recovery_date`.

## Details

Data are generated to match the model's daily interval-censoring: each
case's true onset falls uniformly within its recorded day, and the event
day is the floor of the continuous onset-plus-delay time. So the
recorded day-level delays are exactly a doubly-interval-censored draw,
which makes the simulator suitable for checking calibration, not only
rough recovery.

## Examples

``` r
simulate_linelist(
  n = 5, cfr = 0.6,
  delay = LogNormal(mean = 12.75, sd = 7)
)
#>   onset_date death_date
#> 1 2026-01-02       <NA>
#> 2 2026-01-09 2026-01-23
#> 3 2026-02-26 2026-03-03
#> 4 2026-02-12 2026-03-01
#> 5 2026-02-02 2026-02-14
```
