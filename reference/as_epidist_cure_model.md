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
