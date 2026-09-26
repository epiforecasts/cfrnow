#' Build model inputs from a line list
#'
#' Classifies each case at an observation cut-off into an observed death (with
#' an interval-censored onset-to-death delay), a resolved non-death, or a
#' right-censored survivor still unresolved at the cut-off, and returns the
#' pieces [fit_cfr()] passes to Stan. A case counts as resolved non-death if it
#' has a `recovery_date` on or before the cut-off (or, in a retrospective fit,
#' if it simply never died); such a case contributes the cure term rather than
#' being censored, so recording recoveries tightens the estimate.
#'
#' Onset is taken over the day-window `[onset_lower, onset_upper]` when those
#' columns are present (defaulting to a one-day window at `onset_date`). Deaths
#' and recoveries are recorded to the day.
#'
#' Records that cannot be used are dropped with a warning: a missing onset, an
#' inverted onset window (`onset_upper < onset_lower`), a death with an
#' impossible onset-to-death delay (negative, or longer than `max_delay`), or a
#' recovery dated before onset. In real time, cases whose onset falls after
#' `obs_time` are not yet known and are excluded with a message.
#'
#' @param last_contact_date Optional column name in `linelist` holding the last
#'   date a case was known to be alive and unresolved (a transfer, a discharge
#'   against advice, the last ward note). A case with no recorded outcome is
#'   censored there instead of at `obs_time`, so the follow-up the model sees
#'   stops where the data do. `NA` means followed to the cut-off, and the column
#'   is ignored for a case whose death or recovery was recorded, which was
#'   followed until that happened.
#' @param linelist A data frame with an `onset_date` column, an optional
#'   `onset_lower`/`onset_upper` onset window, a `death_date` column (`NA` for
#'   cases that have not died; use the date the death was notified, i.e. when it
#'   entered the data, so real-time censoring absorbs any reporting lag), and an
#'   optional
#'   `recovery_date` column (`NA` unless the case is a recorded non-fatal
#'   recovery). Dates may be `Date` or coercible.
#' @param obs_time Real-time cut-off: one `Date` (or coercible), one per row of
#'   `linelist`, or the name of a `linelist` column holding them, for data that
#'   reaches the analyst at different times by site. `NULL` gives a
#'   retrospective fit in which every recorded death counts and survivors are
#'   treated as fully resolved. In real time, a case with a recovery on or
#'   before `obs_time` is resolved; one still alive and unresolved is
#'   right-censored; and a death dated after `obs_time` is treated as
#'   not-yet-known (right-censored).
#' @param covariates Character vector of `linelist` column names to carry
#'   through to the per-case model rows, so they can be used in a `prob ~ ...`
#'   formula. The onset date is always carried as `onset`; for a time-varying
#'   CFR derive a time term from it (e.g. `week`) and pass `prob ~ s(week)`.
#' @param t0 Optional time origin (`Date`). Defaults to
#'   `min(onset) - max_delay`.
#' @param max_delay Plausibility filter for data-entry errors, in days: a death
#'   record implying a negative onset-to-death delay, or one longer than
#'   `max_delay`, is dropped as a likely mis-keyed date. This only screens
#'   records; it does **not** bound or truncate the onset-to-death delay the
#'   model fits, so set it comfortably above the longest credible delay to
#'   avoid discarding genuine long-delay deaths (which would bias the delay
#'   short). It also sets the default origin, `t0 = min(onset) - max_delay`.
#'
#' @return A `cfrnow_data` list with the aggregated model inputs (`n_death`,
#'   `death_delay`, `death_width`, `n_recovery`, `recovery_delay`,
#'   `recovery_width`, `n_cens`, `censor_time`, `censor_width`, `n_resolved`,
#'   `n_cases`, `n_deaths`, `n_recoveries`, `t0`, `obs_time`) and a `cases`
#'   data frame with one row per kept case (`y`, `outcome`, `pwindow`,
#'   `swindow`, `onset`, that case's `obs_time`, `follow_up` (days watched,
#'   `Inf` in a retrospective fit) and any requested `covariates`), which
#'   [as_epidist_cure_model()] turns into the model frame.
#' @examples
#' ll <- simulate_linelist(n = 50, delay = LogNormal(2.4, 0.5))
#' prepare_cfr_data(ll, obs_time = as.Date("2026-02-01"))
#' @export
prepare_cfr_data <- function(linelist, obs_time = NULL,
                             covariates = character(), t0 = NULL,
                             max_delay = 60, last_contact_date = NULL) {
  if (!"onset_date" %in% names(linelist)) {
    stop("`linelist` needs an `onset_date` column.", call. = FALSE)
  }
  if (!"death_date" %in% names(linelist)) {
    stop("`linelist` needs a `death_date` column (NA for non-fatal cases).",
      call. = FALSE
    )
  }

  optional_date_col <- function(col, default) {
    if (col %in% names(linelist)) as.Date(linelist[[col]]) else default
  }
  onset <- as.Date(linelist$onset_date)
  onset_lo <- optional_date_col("onset_lower", onset)
  onset_up <- optional_date_col("onset_upper", onset)
  death <- as.Date(linelist$death_date)
  no_recovery <- as.Date(rep(NA, nrow(linelist)))
  recovery <- optional_date_col("recovery_date", no_recovery)
  last_contact <- if (is.null(last_contact_date)) {
    no_recovery
  } else {
    if (!last_contact_date %in% names(linelist)) {
      stop("`last_contact_date` (", last_contact_date, ") is not a column of ",
        "`linelist`.",
        call. = FALSE
      )
    }
    as.Date(linelist[[last_contact_date]])
  }

  # NULL means a retrospective fit. Store the cut-off as a typed Date so
  # downstream code (and the fit) always sees a Date; reject a stray NA so it is
  # not mistaken for a real cut-off. One cut-off per case is allowed, for data
  # that reaches the analyst site by site.
  retrospective <- is.null(obs_time)
  obs_time <- .case_obs_time(obs_time, linelist)
  if (is.null(t0)) t0 <- min(onset, na.rm = TRUE) - max_delay
  t0 <- as.Date(t0)

  onset_lo_day <- as.numeric(onset_lo - t0)
  width <- as.numeric(onset_up - onset_lo) + 1 # onset-window width, days
  # one offset per case: each case is observed up to its own cut-off
  obs_offset <- if (retrospective) {
    rep(Inf, nrow(linelist))
  } else {
    as.numeric(obs_time - t0)
  }

  death_day <- as.numeric(death - t0) # NA for non-fatal
  is_death <- !is.na(death_day) & death_day <= obs_offset
  delay <- death_day - onset_lo_day # from onset-window start
  # Recovered by the cut-off (and not a death by the cut-off) = resolved.
  recovery_day <- as.numeric(recovery - t0)
  recovered <- !is.na(recovery_day) & recovery_day <= obs_offset & !is_death

  # Unusable / erroneous records: missing onset, inverted onset window, a death
  # with an impossible onset->death delay, or a recovery before onset. Guards
  # keep `bad` free of NA.
  bad <- is.na(onset_lo_day) | is.na(width) | width < 1 |
    (is_death & !is.na(delay) & (delay < 0 | delay > max_delay)) |
    (recovered & recovery_day < onset_lo_day)
  n_dropped <- sum(bad)
  if (n_dropped > 0) {
    warning(n_dropped, " unusable record(s) dropped ",
      "(missing onset, inverted window, bad onset-to-death delay, ",
      "or recovery before onset)",
      call. = FALSE
    )
  }

  # Real-time: a case whose onset window opens after the cut-off is not yet
  # known, so exclude it rather than let it inflate the case count.
  future <- !retrospective & !is.na(onset_lo_day) & onset_lo_day > obs_offset
  n_future <- sum(future & !bad)
  if (n_future > 0) {
    message(n_future, " case(s) with onset after the cut-off excluded")
  }

  # Follow-up runs to the end of the cut-off day, or to the end of the last day
  # the case was known unresolved when that comes first. A case whose outcome
  # was recorded was followed until it happened, whatever a last-contact column
  # says: such a column usually carries the outcome's own date, which would
  # otherwise cut the case's follow-up back to the delay it just reported.
  last_contact_day <- as.numeric(last_contact - t0)
  follow_up_end <- obs_offset + 1
  has_contact <- !is.na(last_contact_day) & !is_death & !recovered
  follow_up_end[has_contact] <- pmin(
    follow_up_end[has_contact], last_contact_day[has_contact] + 1
  )

  keep <- !bad & !future
  is_death <- is_death & keep
  recovered <- recovered & keep
  is_surv <- keep & !is_death # every kept non-death

  death_delay <- as.integer(round(delay[is_death]))
  death_width <- width[is_death]

  # Recovered-by-cut-off cases carry their onset-to-recovery timing.
  recovery_delay <- as.integer(round((recovery_day - onset_lo_day)[recovered]))
  recovery_width <- width[recovered]

  # A case with a last-contact date stops being observed there, so it is
  # censored even in a retrospective fit, where the other non-deaths count as
  # resolved.
  lost <- has_contact & is_surv & !recovered
  if (retrospective) {
    # Non-deaths with no recorded recovery are resolved but untimed.
    n_resolved <- sum(is_surv & !recovered & !lost)
    censor_time <- pmax(follow_up_end[lost] - onset_lo_day[lost], 0)
    censor_width <- width[lost]
  } else {
    # Non-deaths not yet recovered are right-censored. A death (or recovery)
    # dated on day `obs_offset` still counts, so the observation horizon is the
    # end of that day, `obs_offset + 1`; a survivor's follow-up runs to there.
    n_resolved <- 0L
    cens <- is_surv & !recovered
    censor_time <- pmax(follow_up_end[cens] - onset_lo_day[cens], 0)
    censor_width <- width[cens]
  }

  # One model row per kept case, carrying the onset date and any requested
  # covariate columns so covariate and time-varying CFR formulas can use them.
  # y is the delay for deaths/recoveries, 0 for resolved non-deaths, and the
  # follow-up time for censored survivors; the likelihood treats each row
  # independently, so row order does not matter.
  n <- nrow(linelist)
  y_case <- integer(n)
  outcome_case <- integer(n)
  outcome_case[is_death] <- .CURE_DEATH
  y_case[is_death] <- as.integer(round(delay[is_death]))
  outcome_case[recovered] <- .CURE_RECOVERY
  y_case[recovered] <-
    as.integer(round((recovery_day - onset_lo_day)[recovered]))
  surv_untimed <- is_surv & !recovered
  censored_case <- if (retrospective) lost else surv_untimed
  resolved_case <- surv_untimed & !censored_case
  outcome_case[resolved_case] <- .CURE_RESOLVED
  outcome_case[censored_case] <- .CURE_CENSORED
  y_case[censored_case] <- as.integer(pmax(
    follow_up_end[censored_case] - onset_lo_day[censored_case], 0
  ))
  # Days of follow-up per kept case, from the start of its onset window to
  # whichever of the cut-off and its last contact comes first. The model reads
  # this as the censored case's `y`; posterior-predictive checks replay it for
  # every case, so a case watched for less time replicates for less time.
  cases <- data.frame(
    y = y_case[keep], outcome = outcome_case[keep],
    pwindow = width[keep], swindow = rep_len(1L, sum(keep)),
    onset = onset[keep], obs_time = obs_time[keep],
    follow_up = pmax(follow_up_end[keep] - onset_lo_day[keep], 0)
  )
  for (cov in covariates) cases[[cov]] <- linelist[[cov]][keep]

  structure(
    list(
      n_death = length(death_delay),
      death_delay = death_delay,
      death_width = death_width,
      n_recovery = length(recovery_delay),
      recovery_delay = recovery_delay,
      recovery_width = recovery_width,
      n_cens = length(censor_time),
      censor_time = censor_time,
      censor_width = censor_width,
      n_resolved = n_resolved,
      n_cases = sum(keep),
      n_deaths = length(death_delay),
      n_recoveries = length(recovery_delay),
      cases = cases,
      t0 = t0,
      # one per kept case, like every other per-case field here; a single
      # cut-off for the whole line list stays a single date
      obs_time = if (length(unique(obs_time[keep])) > 1) {
        obs_time[keep]
      } else {
        obs_time[keep][1]
      }
    ),
    class = "cfrnow_data"
  )
}

#' Per-case observation cut-off
#'
#' `obs_time` is a single date, one date per row of `linelist`, or the name of
#' a `linelist` column holding them; `NULL` is a retrospective fit, which has no
#' cut-off. Returns one date per row either way.
#' @param obs_time The `obs_time` argument as supplied.
#' @param linelist The line list, for a column name and its row count.
#' @return A `Date` vector with one element per row of `linelist`.
#' @noRd
.case_obs_time <- function(obs_time, linelist) {
  n <- nrow(linelist)
  if (is.null(obs_time)) {
    return(as.Date(rep(NA, n)))
  }
  if (is.character(obs_time) && length(obs_time) == 1 &&
        obs_time %in% names(linelist)) {
    obs_time <- linelist[[obs_time]]
  }
  obs_time <- as.Date(obs_time)
  if (!length(obs_time) %in% c(1, n)) {
    stop("`obs_time` must be one date, one per row of `linelist`, or the name ",
      "of a `linelist` column holding them.",
      call. = FALSE
    )
  }
  if (anyNA(obs_time)) {
    stop("`obs_time` must be a valid date, or NULL for a retrospective fit.",
      call. = FALSE
    )
  }
  rep_len(obs_time, n)
}
