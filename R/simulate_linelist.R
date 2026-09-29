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

  otd <- sample_delay(n, delay)
  death_date <- as.Date(rep(NA, n))
  death_date[fatal] <- onset[fatal] + otd[fatal]
  out <- data.frame(onset_date = onset, death_date = death_date)

  if (!is.null(recovery)) {
    otr <- sample_delay(n, recovery)
    recovery_date <- as.Date(rep(NA, n))
    recovery_date[!fatal] <- onset[!fatal] + otr[!fatal]
    out$recovery_date <- recovery_date
  }
  out
}

#' Draw recorded delays from a distspec distribution with fixed parameters
#'
#' Used by [simulate_linelist()] for the onset-to-death and onset-to-recovery
#' delays. The true onset is uniform within its recorded day, so a recorded
#' delay is `floor(onset_frac + delay)` whole days, a doubly-interval-censored
#' draw. A delay with a `max` bounds that recorded delay, as the likelihood
#' [fit_cfr()] uses does, so the offset and the delay are drawn together and
#' resampled as a pair until they fall inside the bound. Errors if any parameter
#' is a prior rather than a fixed number.
#' @param n Number of delays to draw.
#' @param delay A distspec delay distribution with fixed parameters.
#' @return An integer-valued vector of `n` recorded delays (whole days).
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
  d <- switch(fam,
    lognormal = list(
      q = stats::qlnorm, p = stats::plnorm,
      a = pars[["meanlog"]], b = pars[["sdlog"]]
    ),
    gamma = list(
      q = stats::qgamma, p = stats::pgamma,
      a = pars[["shape"]], b = pars[["rate"]]
    ),
    weibull = list(
      q = stats::qweibull, p = stats::pweibull,
      a = pars[["shape"]], b = pars[["scale"]]
    )
  )
  if (is.null(d)) {
    stop("simulate_linelist() supports LogNormal(), Gamma() and Weibull() ",
      "delays only.",
      call. = FALSE
    )
  }
  .draw_recorded(n, d, .delay_max(delay))
}

#' Draw recorded whole-day delays under a bound
#'
#' A recorded delay is `floor(onset_frac + delay)`, and a bound applies to that
#' sum, so offset and delay have to be drawn from their joint distribution
#' conditioned on staying inside it. Conditioning the delay alone would leave
#' the offset uniform, where the conditioning tilts it towards the early part
#' of the day by a factor proportional to `F(max - offset)`.
#'
#' The offset is therefore drawn from that tilted density, by rejection against
#' its largest value, `F(max)`; the delay then follows by inverting its CDF
#' over `[0, F(max - offset)]`. Acceptance depends only on how much of the
#' delay's mass sits in the last day before the bound, never on how far the
#' distribution runs past it, so a delay whose bulk is well beyond the bound
#' costs no more than one whose bulk is inside.
#'
#' @param n Number of delays to draw.
#' @param d A list of the family's quantile and distribution functions (`q`,
#'   `p`) and its two parameters (`a`, `b`), each a scalar or a vector of `n`.
#' @param delay_max The bound, or `Inf`.
#' @return An integer-valued vector of `n` recorded delays (whole days).
#' @noRd
.draw_recorded <- function(n, d, delay_max) {
  frac <- stats::runif(n)
  if (is.infinite(delay_max)) {
    return(floor(frac + d$q(stats::runif(n), d$a, d$b)))
  }
  # a and b are one value, or one per draw
  at <- function(v, i) if (length(v) == 1) v else v[i]
  ceiling_mass <- d$p(delay_max, d$a, d$b)
  if (any(ceiling_mass <= 0)) {
    stop("the delay's `max` (", delay_max, " days) leaves it no probability; ",
      "raise the max, or widen the delay.",
      call. = FALSE
    )
  }
  # Only the offsets still to be accepted are redrawn; testing the settled ones
  # again would need every draw to accept at once, which never happens.
  pending <- seq_len(n)
  while (length(pending) > 0) {
    frac[pending] <- stats::runif(length(pending))
    accept <- stats::runif(length(pending)) * at(ceiling_mass, pending) <=
      d$p(delay_max - frac[pending], at(d$a, pending), at(d$b, pending))
    pending <- pending[!accept]
  }
  u <- stats::runif(n) * d$p(delay_max - frac, d$a, d$b)
  floor(frac + d$q(u, d$a, d$b))
}
