# Build an `epidist_cure_model` from a line list or `prepare_cfr_data()` output

Produces the one-row-per-case data frame the model needs: `y` (a delay
for deaths and recoveries, the follow-up time for survivors), an
`outcome` code, and the primary (`pwindow`) and secondary (`swindow`)
censoring widths. When the input carries timed recoveries the result is
flagged for the two-outcome fit via a `use_recovery` attribute.

## Usage

``` r
as_epidist_cure_model(data)
```

## Arguments

- data:

  Either the list returned by
  [`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
  or a data frame already carrying `y`, `outcome`, `pwindow` and
  `swindow`.

## Value

A data frame of class `epidist_cure_model`.

## Examples

``` r
ll <- simulate_linelist(n = 100, cfr = 0.4, delay = LogNormal(2.4, 0.5))
d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#> 3 case(s) with onset after the cut-off excluded
cure <- as_epidist_cure_model(d)
head(cure)
#>    y outcome pwindow swindow      onset   obs_time follow_up
#> 1 11       0       1       1 2026-02-14 2026-02-24        11
#> 2 33       0       1       1 2026-01-23 2026-02-24        33
#> 3 44       0       1       1 2026-01-12 2026-02-24        44
#> 4 11       1       1       1 2026-02-06 2026-02-24        19
#> 5 18       0       1       1 2026-02-07 2026-02-24        18
#> 6  9       0       1       1 2026-02-16 2026-02-24         9
table(cure$outcome)
#> 
#>  0  1 
#> 68 29 
```
