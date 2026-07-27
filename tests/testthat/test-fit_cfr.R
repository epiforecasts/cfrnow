skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()), error = function(e) FALSE)
  testthat::skip_if(!isTRUE(ok), "cmdstan not installed")
}

# shared sampler settings for the (slow) end-to-end fits
fit_quick <- function(d, ...) {
  fit_cfr(d,
    backend = "cmdstanr", chains = 2, iter = 800, warmup = 400,
    refresh = 0, seed = 1, ...
  )
}

otd <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))

test_that("a retrospective fit matches the naive proportion", {
  skip_if_no_cmdstan()
  set.seed(11)
  ll <- simulate_linelist(n = 1500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL) # every case resolved
  s <- summary(fit_quick(d, delay = otd, cfr_prior = Beta(1, 1)))
  # with every case resolved the cure model reduces to deaths / cases
  expect_equal(s$q50[s$quantity == "cfr"], attr(s, "naive_cfr"), tolerance = 0.03)
})

test_that("real-time correction lifts above naive and covers truth", {
  skip_if_no_cmdstan()
  set.seed(12)
  ll <- simulate_linelist(
    n = 2000, cfr = 0.6, onset_days = 40,
    delay = LogNormal(2.4, 0.5)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  s <- summary(fit_quick(d, delay = otd, cfr_prior = Beta(1, 1)))
  cfr <- s[s$quantity == "cfr", ]
  expect_gt(cfr[["q50"]], attr(s, "naive_cfr")) # corrected > naive
  expect_lt(cfr[["q2.5"]], 0.6) # 95% CrI covers the truth
  expect_gt(cfr[["q97.5"]], 0.6)
  expect_true(all(s$rhat < 1.05))
})

test_that("a gamma delay recovers the CFR and delay moments", {
  skip_if_no_cmdstan()
  set.seed(2)
  ll <- simulate_linelist(
    n = 2000, cfr = 0.4, onset_days = 40,
    delay = Gamma(mean = 8, sd = 4)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  s <- summary(fit_quick(d,
    delay = Gamma(
      shape = Normal(4, 1),
      rate = Normal(0.5, 0.2)
    ),
    cfr_prior = Beta(1, 1)
  ))
  expect_lt(s[s$quantity == "cfr", "q2.5"], 0.4)
  expect_gt(s[s$quantity == "cfr", "q97.5"], 0.4)
  expect_equal(s[s$quantity == "delay_mean", "q50"], 8, tolerance = 0.8)
})

test_that("a weibull delay recovers the CFR and delay moments", {
  skip_if_no_cmdstan()
  set.seed(3)
  ll <- simulate_linelist(
    n = 2000, cfr = 0.4, onset_days = 40,
    delay = Weibull(shape = 1.5, scale = 13)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  s <- summary(fit_quick(d,
    delay = Weibull(shape = Normal(1.5, 0.4), scale = Normal(13, 3)),
    cfr_prior = Beta(1, 1)
  ))
  true_mean <- 13 * gamma(1 + 1 / 1.5) # Weibull mean = scale * Gamma(1 + 1/shape)
  expect_lt(s[s$quantity == "cfr", "q2.5"], 0.4)
  expect_gt(s[s$quantity == "cfr", "q97.5"], 0.4)
  expect_equal(s[s$quantity == "delay_mean", "q50"], true_mean, tolerance = 1.5)
})

test_that("a Weibull recovery delay compiles and gives sensible moments", {
  skip_if_no_cmdstan()
  set.seed(5)
  ll <- simulate_linelist(
    n = 2000, cfr = 0.4, onset_days = 40,
    delay = Gamma(mean = 12.75, sd = 7),
    recovery = Weibull(mean = 21, sd = 9)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  fit <- fit_quick(d,
    delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)),
    recovery_delay = Weibull(shape = Normal(2.4, 0.6), scale = Normal(23.7, 4)),
    cfr_prior = Beta(1, 1)
  )
  expect_equal(fit$cfrnow$recovery_family, "weibull") # differs from gamma death
  s <- summary(fit)
  rm <- s[s$quantity == "recovery_mean", ] # onset-to-recovery ~21 (heavy censoring)
  expect_gt(rm[["q50"]], 17)
  expect_lt(rm[["q50"]], 26)
})

test_that("a fixed delay runs the Ghani/Nishiura estimator (delay held constant)", {
  skip_if_no_cmdstan()
  set.seed(4)
  ll <- simulate_linelist(n = 1500, cfr = 0.5, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  s <- summary(fit_quick(d,
    delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
    cfr_prior = Beta(1, 1)
  ))
  expect_true(is.na(s[s$quantity == "delay_mean", "rhat"])) # delay is constant
  expect_lt(s[s$quantity == "cfr", "q2.5"], 0.5)
  expect_gt(s[s$quantity == "cfr", "q97.5"], 0.5)
})

test_that("the formula interface puts covariates on the CFR", {
  skip_if_no_cmdstan()
  set.seed(20)
  a <- simulate_linelist(n = 1200, cfr = 0.2, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 1200, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  da <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  db <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  da$grp <- "low"
  db$grp <- "high"
  d <- as_epidist_cure_model(rbind(da, db))
  fit <- fit_quick(d, delay = otd, formula = brms::bf(mu ~ 1, cfr ~ grp))
  fe <- brms::fixef(fit)
  expect_true("cfr_grplow" %in% rownames(fe))
  expect_lt(fe["cfr_grplow", "Estimate"], 0) # low group < high group
})

test_that("a two-outcome fit uses recovery timing (own family) and recovers F_R", {
  skip_if_no_cmdstan()
  set.seed(5)
  ll <- simulate_linelist(
    n = 2000, cfr = 0.4, onset_days = 40,
    delay = Gamma(mean = 12.75, sd = 7),
    recovery = LogNormal(mean = 21, sd = 9)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  fit <- fit_quick(d,
    delay = Gamma(shape = Normal(3.3, 1), rate = Normal(0.26, 0.08)),
    recovery_delay = LogNormal(meanlog = Normal(2.9, 0.3), sdlog = Normal(0.5, 0.2)),
    cfr_prior = Beta(1, 1)
  )
  expect_equal(fit$cfrnow$recovery_family, "lognormal") # differs from gamma death
  s <- summary(fit)
  rm <- s[s$quantity == "recovery_mean", ] # true 21
  expect_lt(rm[["q2.5"]], 21)
  expect_gt(rm[["q97.5"]], 21)
})

test_that("a young outbreak is flagged low-information and print warns", {
  skip_if_no_cmdstan()
  set.seed(10)
  ll <- simulate_linelist(
    n = 25, cfr = 0.5, onset_days = 10,
    delay = LogNormal(2.4, 0.5)
  )
  d <- prepare_cfr_data(ll, obs_time = min(ll$onset_date) + 3)
  fit <- fit_quick(d,
    delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
    cfr_prior = Beta(1, 1)
  )
  expect_true(attr(summary(fit), "cfr_low_information"))
  expect_message(print(fit), "weakly identified")
})

test_that("print reports the delay family, counts and naive CFR", {
  skip_if_no_cmdstan()
  set.seed(6)
  ll <- simulate_linelist(n = 400, cfr = 0.5, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_quick(d, delay = otd, cfr_prior = Beta(1, 1))
  expect_message(print(fit), "lognormal delay")
  expect_message(print(fit), "naive CFR")
  expect_invisible(suppressMessages(print(fit)))
})

test_that("an intercept-free cfr formula fits one CFR per group", {
  skip_if_no_cmdstan()
  set.seed(21)
  a <- simulate_linelist(n = 1200, cfr = 0.25, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 1200, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  da <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  db <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  da$grp <- "low"
  db$grp <- "high"
  d <- as_epidist_cure_model(rbind(da, db))
  fit <- fit_quick(d, delay = otd, formula = brms::bf(mu ~ 1, cfr ~ 0 + grp))
  fe <- brms::fixef(fit)
  expect_true(all(c("cfr_grplow", "cfr_grphigh") %in% rownames(fe)))
  expect_lt(fe["cfr_grplow", "Estimate"], fe["cfr_grphigh", "Estimate"])
})

test_that("summary() reports a CFR per group for a grouped fit", {
  skip_if_no_cmdstan()
  set.seed(22)
  a <- simulate_linelist(n = 800, cfr = 0.25, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 800, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  da <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  db <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  da$site <- "A"
  db$site <- "B"
  d <- as_epidist_cure_model(rbind(da, db))
  s <- summary(fit_quick(d, delay = otd, formula = brms::bf(mu ~ 1, cfr ~ 0 + site)))
  expect_true(all(c("cfr[A]", "cfr[B]") %in% s$quantity))
  expect_lt(abs(s[s$quantity == "cfr[A]", "q50"] - 0.25), 0.07)
  expect_lt(abs(s[s$quantity == "cfr[B]", "q50"] - 0.60), 0.07)
  # the weak-identification flag is not defined for a grouped fit
  expect_true(is.na(attr(s, "cfr_low_information")))
})
