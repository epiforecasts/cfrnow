# Build model inputs from a line list

Classifies each case at an observation cut-off into an observed death
(with an interval-censored onset-to-death delay), a resolved non-death,
or a right-censored survivor still unresolved at the cut-off, and
returns the pieces
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md)
passes to Stan. A case counts as resolved non-death if it has a
`recovery_date` on or before the cut-off (or, in a retrospective fit, if
it simply never died); such a case contributes the cure term rather than
being censored, so recording recoveries tightens the estimate.

## Usage

``` r
prepare_cfr_data(
  linelist,
  obs_time = NULL,
  covariates = character(),
  t0 = NULL,
  max_delay = 60
)
```

## Arguments

- linelist:

  A data frame with an `onset_date` column, an optional
  `onset_lower`/`onset_upper` onset window, a `death_date` column (`NA`
  for cases that have not died; use the date the death was notified,
  i.e. when it entered the data, so real-time censoring absorbs any
  reporting lag), and an optional `recovery_date` column (`NA` unless
  the case is a recorded non-fatal recovery). Dates may be `Date` or
  coercible.

- obs_time:

  Real-time cut-off (`Date` or coercible), or `NULL` for a retrospective
  fit in which every recorded death counts and survivors are treated as
  fully resolved. In real time, a case with a recovery on or before
  `obs_time` is resolved; one still alive and unresolved is
  right-censored; and a death dated after `obs_time` is treated as
  not-yet-known (right-censored).

- covariates:

  Character vector of `linelist` column names to carry through to the
  per-case model rows, so they can be used in a `cfr ~ ...` formula. The
  onset date is always carried as `onset`; for a time-varying CFR derive
  a time term from it (e.g. `week`) and pass `cfr ~ s(week)`.

- t0:

  Optional time origin (`Date`). Defaults to `min(onset) - max_delay`.

- max_delay:

  Plausibility filter for data-entry errors, in days: a death record
  implying a negative onset-to-death delay, or one longer than
  `max_delay`, is dropped as a likely mis-keyed date. This only screens
  records; it does **not** bound or truncate the onset-to-death delay
  the model fits, so set it comfortably above the longest credible delay
  to avoid discarding genuine long-delay deaths (which would bias the
  delay short). It also sets the default origin,
  `t0 = min(onset) - max_delay`.

## Value

A `cfrnow_data` list with the aggregated model inputs (`n_death`,
`death_delay`, `death_width`, `n_recovery`, `recovery_delay`,
`recovery_width`, `n_cens`, `censor_time`, `censor_width`, `n_resolved`,
`n_cases`, `n_deaths`, `n_recoveries`, `t0`, `obs_time`) and a `cases`
data frame with one row per kept case (`y`, `outcome`, `pwindow`,
`swindow`, `onset` and any requested `covariates`), which
[`as_epidist_cure_model()`](https://epiforecasts.io/cfrnow/reference/as_epidist_cure_model.md)
turns into the model frame.

## Details

Onset is taken over the day-window `[onset_lower, onset_upper]` when
those columns are present (defaulting to a one-day window at
`onset_date`). Deaths and recoveries are recorded to the day.

Records that cannot be used are dropped with a warning: a missing onset,
an inverted onset window (`onset_upper < onset_lower`), a death with an
impossible onset-to-death delay (negative, or longer than `max_delay`),
or a recovery dated before onset. In real time, cases whose onset falls
after `obs_time` are not yet known and are excluded with a message.

## Examples

``` r
ll <- simulate_linelist(n = 50, delay = LogNormal(2.4, 0.5))
prepare_cfr_data(ll, obs_time = as.Date("2026-02-01"))
#> 21 case(s) with onset after the cut-off excluded
#> $n_death
#> [1] 8
#> 
#> $death_delay
#> [1] 14 19 11 14 26  9  8  5
#> 
#> $death_width
#> [1] 1 1 1 1 1 1 1 1
#> 
#> $n_recovery
#> [1] 0
#> 
#> $recovery_delay
#> integer(0)
#> 
#> $recovery_width
#> numeric(0)
#> 
#> $n_cens
#> [1] 21
#> 
#> $censor_time
#>  [1] 10 21  2 29 24  9 28  8 11 26  1 27 24 11 23  1 16 11  2 32 27
#> 
#> $censor_width
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $n_resolved
#> [1] 0
#> 
#> $n_cases
#> [1] 29
#> 
#> $n_deaths
#> [1] 8
#> 
#> $n_recoveries
#> [1] 0
#> 
#> $cases
#>     y outcome pwindow swindow      onset
#> 1  10       0       1       1 2026-01-23
#> 2  21       0       1       1 2026-01-12
#> 3   2       0       1       1 2026-01-31
#> 4  29       0       1       1 2026-01-04
#> 5  24       0       1       1 2026-01-09
#> 6  14       1       1       1 2026-01-05
#> 7  19       1       1       1 2026-01-05
#> 8   9       0       1       1 2026-01-24
#> 9  11       1       1       1 2026-01-15
#> 10 14       1       1       1 2026-01-13
#> 11 26       1       1       1 2026-01-02
#> 12 28       0       1       1 2026-01-05
#> 13  9       1       1       1 2026-01-21
#> 14  8       1       1       1 2026-01-06
#> 15  5       1       1       1 2026-01-04
#> 16  8       0       1       1 2026-01-25
#> 17 11       0       1       1 2026-01-22
#> 18 26       0       1       1 2026-01-07
#> 19  1       0       1       1 2026-02-01
#> 20 27       0       1       1 2026-01-06
#> 21 24       0       1       1 2026-01-09
#> 22 11       0       1       1 2026-01-22
#> 23 23       0       1       1 2026-01-10
#> 24  1       0       1       1 2026-02-01
#> 25 16       0       1       1 2026-01-17
#> 26 11       0       1       1 2026-01-22
#> 27  2       0       1       1 2026-01-31
#> 28 32       0       1       1 2026-01-01
#> 29 27       0       1       1 2026-01-06
#> 
#> $t0
#> [1] "2025-11-02"
#> 
#> $obs_time
#> [1] "2026-02-01"
#> 
#> attr(,"class")
#> [1] "cfrnow_data"
```
