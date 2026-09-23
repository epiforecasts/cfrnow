# Precompute the time-varying CFR demonstration used in the vignette.
#
# The vignette cannot run Stan at build time (no CmdStan on the builder), so we
# fit here, save the recovered cfr(week) posterior summary alongside the known
# true curve, and let the vignette read the cached result and plot it.
#
# Run from the package root with CmdStan available:
#   Rscript data-raw/time_varying_cfr.R

library(cfrnow)

set.seed(20260709)

# --- A line list with a genuinely time-varying true CFR -------------------
# Onsets span twelve weeks; the true CFR ramps logistically from ~0.15 to ~0.6
# over that window, so a smooth on the onset week has a real signal to recover.
n <- 4000
onset_start <- as.Date("2026-01-01")
onset_days <- 12 * 7

true_cfr <- function(week) stats::plogis(-1.7 + 0.19 * week)

onset_day <- sample.int(onset_days, n, replace = TRUE) - 1
onset_date <- onset_start + onset_day
onset_week <- onset_day %/% 7

fatal <- stats::runif(n) < true_cfr(onset_week)

# Onset-to-death delay, interval-censored to the day exactly as the model
# assumes (true onset uniform within its recorded day; event day is the floor).
onset_frac <- stats::runif(n)
otd <- stats::rlnorm(n, meanlog = 2.41, sdlog = 0.51)
death_date <- as.Date(rep(NA, n))
death_date[fatal] <- onset_date[fatal] + floor(onset_frac[fatal] + otd[fatal])

ll <- data.frame(onset_date = onset_date, death_date = death_date)

# --- Prepare at a real-time cut-off and derive the week index -------------
# The cut-off sits three weeks past the last onset: older weeks are mostly
# resolved, the most recent weeks still carry censored survivors, so the
# recovered curve should widen towards the present.
cutoff <- max(ll$onset_date) + 21
cure <- as_epidist_cure_model(prepare_cfr_data(ll, obs_time = cutoff))
cure$week <- as.numeric(cure$onset - min(cure$onset)) %/% 7

# A fixed (Ghani/Nishiura) onset-to-death delay isolates the time-varying-CFR
# demonstration from the delay's uncertainty.
onset_to_death <- LogNormal(meanlog = 2.41, sdlog = 0.51)

# A fixed-df natural spline on the week enters `prob` as ordinary basis
# columns with Normal priors, so there is no penalised-smooth variance
# hyperparameter and hence no funnel: the model samples cleanly at the default
# adapt_delta while still fitting a flexible curve.
fit <- fit_cfr(
  cure,
  delay = onset_to_death,
  prob_prior = Beta(1, 1),
  formula = brms::bf(mu ~ 1, prob ~ splines::ns(week, df = 3)),
  backend = "cmdstanr", chains = 4, cores = 4, iter = 1000, refresh = 0,
  seed = 1
)

# --- Recover prob(week) on the prob dpar (logit link) ----------------------
# Predict on the fitted rows rather than a fresh grid: every case in a week
# shares the same `week`, so its prob prediction is identical, and reusing the
# training design matrix sidesteps recomputing the spline basis from a new
# range.
weeks <- sort(unique(cure$week))
cfr_all <- brms::posterior_epred(fit, dpar = "prob")   # ndraws x n_cases
cfr_draws <- cfr_all[, match(weeks, cure$week), drop = FALSE]

cfr_week <- data.frame(
  week = weeks,
  n_cases = as.integer(table(factor(cure$week, levels = weeks))),
  true_cfr = true_cfr(weeks),
  median = apply(cfr_draws, 2, stats::median),
  lower = apply(cfr_draws, 2, stats::quantile, probs = 0.05),
  upper = apply(cfr_draws, 2, stats::quantile, probs = 0.95)
)

saveRDS(cfr_week, file.path("inst", "vignette-data", "time_varying_cfr.rds"))
print(cfr_week, row.names = FALSE)
