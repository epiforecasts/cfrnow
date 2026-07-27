#' Simulate a line list for testing and examples
#'
#' Draws onset dates over a window, marks each case fatal with probability
#' `cfr`, and gives fatal cases an onset-to-death delay drawn from `delay`. Pass
#' a `recovery` delay to also simulate onset-to-recovery times for the non-fatal
#' cases and add a `recovery_date` column. Delays are distspec distributions
#' with fixed parameters, matching what [fit_cfr()] takes; you can specify them
#' by mean and sd (e.g. `LogNormal(mean = 12.75, sd = 7)`). Returns a
#' line list with `onset_date`, `death_date` (`NA` for non-fatal cases) and,
#' when `recovery` is given, `recovery_date` (`NA` for fatal cases). The full,
#' untruncated outcomes are simulated; pass the result to [prepare_cfr_data()]
#' with an `obs_time` to induce the real-time truncation.
#'
#' Data are generated to match the model's daily interval-censoring: each case's
#' true onset falls uniformly within its recorded day, and the event day is the
#' floor of the continuous onset-plus-delay time. So the recorded day-level
#' delays are exactly a doubly-interval-censored draw, which makes the simulator
#' suitable for checking calibration, not only rough recovery.
#'
#' @param n Number of cases.
#' @param cfr True case fatality ratio.
#' @param delay Onset-to-death delay: a distspec distribution
#'   ([distspec::LogNormal()], [distspec::Gamma()] or [distspec::Weibull()])
#'   with fixed parameters.
#' @param recovery Optional onset-to-recovery delay (same form as `delay`); when
#'   given, non-fatal cases get a `recovery_date`.
#' @param onset_start First possible onset date.
#' @param onset_days Width of the onset window (days); onsets are uniform
#'   over it.
#' @return A data frame with `onset_date`, `death_date` and, if `recovery` is
#'   given, `recovery_date`.
#' @examples
#' simulate_linelist(
#'   n = 5, cfr = 0.6,
#'   delay = LogNormal(mean = 12.75, sd = 7)
#' )
#' @export
simulate_linelist <- function(n = 200, cfr = 0.5, delay, recovery = NULL,
                              onset_start = as.Date("2026-01-01"),
                              onset_days = 60) {
  if (missing(delay)) {
    stop("supply a `delay` (a distspec distribution with fixed parameters).",
      call. = FALSE
    )
  }
  onset <- as.Date(onset_start) + sample.int(onset_days, n, replace = TRUE) - 1
  fatal <- stats::runif(n) < cfr
  # Sub-day onset offset: true onset is uniform within its recorded day.
  onset_frac <- stats::runif(n)

  otd <- sample_delay(n, delay)
  death_date <- as.Date(rep(NA, n))
  # Event day = floor of the continuous onset-plus-delay time, so day-level
  # delays are a genuine doubly-interval-censored draw.
  death_date[fatal] <- onset[fatal] + floor(onset_frac[fatal] + otd[fatal])
  out <- data.frame(onset_date = onset, death_date = death_date)

  if (!is.null(recovery)) {
    otr <- sample_delay(n, recovery)
    recovery_date <- as.Date(rep(NA, n))
    recovery_date[!fatal] <-
      onset[!fatal] + floor(onset_frac[!fatal] + otr[!fatal])
    out$recovery_date <- recovery_date
  }
  out
}

#' Draw delays from a distspec distribution with fixed parameters
#'
#' Used by [simulate_linelist()] for the onset-to-death and onset-to-recovery
#' delays. Errors if any parameter is a prior rather than a fixed number.
#' @param n Number of delays to draw.
#' @param delay A distspec delay distribution with fixed parameters.
#' @return A numeric vector of `n` delays (days).
#' @noRd
sample_delay <- function(n, delay) {
  fam <- get_distribution(delay)
  pars <- get_parameters(delay)[natural_params(delay)]
  if (!all(vapply(pars, is.numeric, logical(1)))) {
    stop("simulate_linelist() needs a delay with fixed parameters (numbers), ",
      "not priors.",
      call. = FALSE
    )
  }
  out <- switch(fam,
    lognormal = stats::rlnorm(n, pars[["meanlog"]], pars[["sdlog"]]),
    gamma = stats::rgamma(n, shape = pars[["shape"]], rate = pars[["rate"]]),
    weibull = stats::rweibull(
      n,
      shape = pars[["shape"]], scale = pars[["scale"]]
    )
  )
  if (is.null(out)) {
    stop("simulate_linelist() supports LogNormal(), Gamma() and Weibull() ",
      "delays only.",
      call. = FALSE
    )
  }
  out
}
