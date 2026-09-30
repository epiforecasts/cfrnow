# Precompute the fits behind the "Real-time CFR over time" vignette.
#
# The vignette builder has no CmdStan, so we fit here and save per-week
# summaries that the vignette reads and plots.
#
# Run from the package root with CmdStan available:
#   Rscript data-raw/real_time_cfr.R

library(cfrnow)

set.seed(20260930)

# --- A line list whose CFR steps up late in the outbreak -------------------
# Twelve onset weeks; the true CFR is 0.25 up to week 7 and 0.5 from week 8,
# so the change falls in the weeks that are least resolved at an early cut-off.
n <- 2000
onset_start <- as.Date("2026-01-01")
true_cfr <- function(week) ifelse(week >= 8, 0.5, 0.25)

onset_day <- sample.int(12 * 7, n, replace = TRUE) - 1
onset_date <- onset_start + onset_day
week <- onset_day %/% 7
fatal <- stats::runif(n) < true_cfr(week)

# Onset-to-death delay, interval-censored to the day as the model assumes
death_date <- as.Date(rep(NA, n))
death_date[fatal] <- onset_date[fatal] +
  floor(stats::runif(sum(fatal)) + stats::rlnorm(sum(fatal), 2.41, 0.51))
ll <- data.frame(onset_date = onset_date, death_date = death_date)

# The fraction of each week's cases that die, eventually: the quantity the
# nowcast estimates for the cases actually seen, as opposed to `true_cfr`
truth <- data.frame(
  week = 0:11,
  true_cfr = true_cfr(0:11),
  realised = as.numeric(tapply(fatal, factor(week, levels = 0:11), mean))
)

# --- Fit at two cut-offs, with a constant and a time-varying prob ----------
onset_to_death <- LogNormal(
  meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)
)
formulas <- list(
  constant = brms::bf(mu ~ 1, prob ~ 1),
  spline = brms::bf(mu ~ 1, prob ~ splines::ns(week, df = 3))
)
cutoffs <- c(early = 7, late = 42) # days after the last onset

summarise_draws <- function(draws, week, what) {
  weeks <- sort(unique(week))
  per_week <- vapply(weeks, function(w) {
    rowMeans(draws[, week == w, drop = FALSE])
  }, numeric(nrow(draws)))
  data.frame(
    week = weeks, quantity = what,
    median = apply(per_week, 2, stats::median),
    lower = apply(per_week, 2, stats::quantile, probs = 0.05),
    upper = apply(per_week, 2, stats::quantile, probs = 0.95)
  )
}

results <- list()
naive <- list()
for (cut in names(cutoffs)) {
  cure <- as_epidist_cure_model(prepare_cfr_data(ll,
    obs_time = max(ll$onset_date) + cutoffs[[cut]]
  ))
  cure$week <- as.numeric(cure$onset - onset_start) %/% 7

  naive[[cut]] <- data.frame(
    cutoff = cut, week = 0:11,
    naive = as.numeric(tapply(
      cure$outcome == 1, factor(cure$week, levels = 0:11), mean
    ))
  )

  for (model in names(formulas)) {
    fit <- fit_cfr(cure,
      delay = onset_to_death, prob_prior = Beta(1, 1),
      formula = formulas[[model]],
      backend = "cmdstanr", chains = 4, cores = 4, iter = 1000,
      refresh = 0, seed = 1
    )
    # the fit stores each fitted case's onset date in row order
    fit_week <- as.numeric(fit$cfrnow$onset - onset_start) %/% 7
    prob <- brms::posterior_epred(fit, dpar = "prob")
    pi <- posterior::as_draws_matrix(posterior_prob_death(fit))
    # one outcome per case and draw, so the weekly fraction varies with whether
    # each unresolved case dies as well as with the parameters
    fatal_draws <- matrix(stats::rbinom(length(pi), 1, pi), nrow(pi))
    results[[paste(cut, model)]] <- cbind(
      cutoff = cut, model = model,
      rbind(
        summarise_draws(prob, fit_week, "prob"),
        summarise_draws(fatal_draws, fit_week, "nowcast")
      )
    )
  }
}

real_time_cfr <- list(
  estimates = do.call(rbind, results),
  naive = do.call(rbind, naive),
  truth = truth
)
rownames(real_time_cfr$estimates) <- NULL
saveRDS(real_time_cfr, file.path("inst", "vignette-data", "real_time_cfr.rds"))
print(real_time_cfr)
