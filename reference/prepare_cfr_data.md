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
  max_delay = 60,
  last_contact_date = NULL
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

  Real-time cut-off, or `NULL` for a retrospective fit in which every
  recorded death counts and survivors are treated as fully resolved. A
  single `Date` (or coercible) is shared by every case; a vector with
  one entry per row of `linelist` gives each case its own cut-off; or a
  string naming a `linelist` column holding those per-case cut-offs. In
  real time, a case with a recovery on or before its own `obs_time` is
  resolved; one still alive and unresolved is right-censored; and a
  death dated after its own `obs_time` is treated as not-yet-known
  (right-censored).

- covariates:

  Character vector of `linelist` column names to carry through to the
  per-case model rows, so they can be used in a `prob ~ ...` formula.
  The onset date is always carried as `onset`; for a time-varying CFR
  derive a time term from it (e.g. `week`) and pass `prob ~ s(week)`.

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

- last_contact_date:

  Optional column name in `linelist` holding the last date a case was
  known to be alive and unresolved (a transfer, a discharge against
  advice, the last ward note). A case with no recorded outcome is
  censored there instead of at `obs_time`, so the follow-up the model
  sees stops where the data do. `NA` means followed to the cut-off, and
  the column is ignored for a case whose death or recovery was recorded,
  which was followed until that happened.

## Value

A `cfrnow_data` list with the aggregated model inputs (`n_death`,
`death_delay`, `death_width`, `n_recovery`, `recovery_delay`,
`recovery_width`, `n_cens`, `censor_time`, `censor_width`, `n_resolved`,
`n_cases`, `n_deaths`, `n_recoveries`, `t0`, `obs_time`, one per kept
case) and a `cases` data frame with one row per kept case (`y`,
`outcome`, `pwindow`, `swindow`, `onset`, `obs_time`, `follow_up` (days
watched, `Inf` in a retrospective fit) and any requested `covariates`),
which
[`as_epidist_cure_model()`](https://epiforecasts.io/cfrnow/reference/as_epidist_cure_model.md)
turns into the model frame.

## Details

Onset is taken over the day-window `[onset_lower, onset_upper]` when
those columns are present, defaulting to a one-day window at
`onset_date` for a case whose window is missing (`NA`), or for every
case when the columns are absent. In real time, a window that closes
after the case's cut-off is cut back to the cut-off, since the case was
already known then. Deaths and recoveries are recorded to the day.

Records that cannot be used are dropped with a warning: a missing onset,
an inverted onset window (`onset_upper < onset_lower`), a death with an
impossible onset-to-death delay (negative, or longer than `max_delay`),
or a recovery dated before onset. In real time, a case whose onset falls
after its own `obs_time` is not yet known and is excluded with a
message.

Data often reaches the analyst at different times by site, so `obs_time`
may give each case its own cut-off rather than one shared by the whole
line list: a case's deaths, recoveries, follow-up and exclusion are all
judged against its own cut-off, not the latest one in the data.

## Examples

``` r
ll <- simulate_linelist(n = 50, delay = LogNormal(2.4, 0.5))
prepare_cfr_data(ll, obs_time = as.Date("2026-02-01"))
#> 21 case(s) with onset after the cut-off excluded
#> $n_death
#> [1] 9
#> 
#> $death_delay
#> [1] 15 20  8  3 11 19  6  5  4
#> 
#> $death_width
#> [1] 1 1 1 1 1 1 1 1 1
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
#> [1] 20
#> 
#> $censor_time
#>  [1] 10 21  2 29 24 20 28 12  8 11 26  1 27 24 23  1 16  2 32 27
#> 
#> $censor_width
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $n_resolved
#> [1] 0
#> 
#> $n_cases
#> [1] 29
#> 
#> $n_deaths
#> [1] 9
#> 
#> $n_recoveries
#> [1] 0
#> 
#> $cases
#>     y outcome pwindow swindow      onset   obs_time follow_up
#> 1  10       0       1       1 2026-01-23 2026-02-01        10
#> 2  21       0       1       1 2026-01-12 2026-02-01        21
#> 3   2       0       1       1 2026-01-31 2026-02-01         2
#> 4  29       0       1       1 2026-01-04 2026-02-01        29
#> 5  24       0       1       1 2026-01-09 2026-02-01        24
#> 6  15       1       1       1 2026-01-05 2026-02-01        28
#> 7  20       1       1       1 2026-01-05 2026-02-01        28
#> 8   8       1       1       1 2026-01-24 2026-02-01         9
#> 9   3       1       1       1 2026-01-15 2026-02-01        18
#> 10 20       0       1       1 2026-01-13 2026-02-01        20
#> 11 11       1       1       1 2026-01-02 2026-02-01        31
#> 12 28       0       1       1 2026-01-05 2026-02-01        28
#> 13 12       0       1       1 2026-01-21 2026-02-01        12
#> 14 19       1       1       1 2026-01-06 2026-02-01        27
#> 15  6       1       1       1 2026-01-04 2026-02-01        29
#> 16  8       0       1       1 2026-01-25 2026-02-01         8
#> 17 11       0       1       1 2026-01-22 2026-02-01        11
#> 18 26       0       1       1 2026-01-07 2026-02-01        26
#> 19  1       0       1       1 2026-02-01 2026-02-01         1
#> 20 27       0       1       1 2026-01-06 2026-02-01        27
#> 21 24       0       1       1 2026-01-09 2026-02-01        24
#> 22  5       1       1       1 2026-01-22 2026-02-01        11
#> 23 23       0       1       1 2026-01-10 2026-02-01        23
#> 24  1       0       1       1 2026-02-01 2026-02-01         1
#> 25 16       0       1       1 2026-01-17 2026-02-01        16
#> 26  4       1       1       1 2026-01-22 2026-02-01        11
#> 27  2       0       1       1 2026-01-31 2026-02-01         2
#> 28 32       0       1       1 2026-01-01 2026-02-01        32
#> 29 27       0       1       1 2026-01-06 2026-02-01        27
#> 
#> $t0
#> [1] "2025-11-02"
#> 
#> $obs_time
#>  [1] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#>  [6] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#> [11] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#> [16] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#> [21] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#> [26] "2026-02-01" "2026-02-01" "2026-02-01" "2026-02-01"
#> 
#> attr(,"class")
#> [1] "cfrnow_data"
```
