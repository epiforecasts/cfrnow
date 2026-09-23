test_that("pp_check_cfr rejects objects that are not cfrnow fits", {
  expect_error(pp_check_cfr(1), "must come from fit_cfr")
  expect_error(
    pp_check_cfr(structure(list(), class = "brmsfit")),
    "must come from fit_cfr"
  )
})

test_that(".cfr_ppc_stats counts deaths and replays truncation", {
  set.seed(1)
  nd <- 20
  n <- 50
  loc <- matrix(2.4, nd, n) # meanlog
  sc <- matrix(0.5, nd, n) # sdlog
  obs <- list(deaths = 30L, recoveries = NA_integer_, death_delays = c(5, 8, 12))

  # every case fatal, no truncation: each draw yields n observed deaths
  all_fatal <- .cfr_ppc_stats(
    matrix(1, nd, n), loc, sc, rep(Inf, n), "lognormal",
    FALSE, NULL, NULL, NULL, obs
  )
  expect_equal(nrow(all_fatal$counts), nd)
  expect_true(all(all_fatal$counts$n == n))
  expect_equal(all_fatal$observed_counts$n, 30L)
  expect_equal(nrow(all_fatal$observed_delays), 3L)

  # no fatal cases: no observed deaths anywhere
  none_fatal <- .cfr_ppc_stats(
    matrix(0, nd, n), loc, sc, rep(Inf, n), "lognormal",
    FALSE, NULL, NULL, NULL, obs
  )
  expect_true(all(none_fatal$counts$n == 0))
  expect_equal(nrow(none_fatal$delays), 0L)

  # a one-day horizon truncates almost every (long-delay) death away
  truncated <- .cfr_ppc_stats(
    matrix(1, nd, n), loc, sc, rep(1, n), "lognormal",
    FALSE, NULL, NULL, NULL, obs
  )
  expect_lt(mean(truncated$counts$n), n)
})

test_that(".cfr_ppc_stats handles gamma delays and the recovery branch", {
  set.seed(2)
  nd <- 15
  n <- 40
  gamma_all <- .cfr_ppc_stats(
    matrix(1, nd, n), matrix(10, nd, n), matrix(3, nd, n), rep(Inf, n),
    "gamma", FALSE, NULL, NULL, NULL,
    list(deaths = 5L, recoveries = NA_integer_, death_delays = 1:3)
  )
  expect_true(all(gamma_all$counts$n == n))

  two_outcome <- .cfr_ppc_stats(
    matrix(0.5, nd, n), matrix(2.4, nd, n), matrix(0.5, nd, n), rep(Inf, n),
    "lognormal", TRUE,
    matrix(2.0, nd, n), matrix(0.4, nd, n), "lognormal",
    list(deaths = 20L, recoveries = 10L, death_delays = 1:4)
  )
  expect_setequal(unique(two_outcome$counts$outcome), c("deaths", "recoveries"))
  expect_setequal(two_outcome$observed_counts$outcome, c("deaths", "recoveries"))
})

test_that("the ppc plot helpers return ggplots", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  reps <- .cfr_ppc_stats(
    matrix(0.4, 20, 50), matrix(2.4, 20, 50), matrix(0.5, 20, 50),
    rep(30, 50), "lognormal", FALSE, NULL, NULL, NULL,
    list(deaths = 15L, recoveries = NA_integer_, death_delays = c(5, 8, 12, 15))
  )
  expect_s3_class(.cfr_ppc_counts_plot(reps), "ggplot")
  expect_s3_class(.cfr_ppc_delay_plot(reps), "ggplot")
})

test_that("pp_check_cfr returns ggplots and reproduces the death count", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(1)
  ll <- simulate_linelist(n = 300, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, seed = 1
  )

  # the cut-off is recorded so the check can replay the real-time truncation
  expect_false(is.na(fit$cfrnow$obs_time))
  expect_s3_class(pp_check_cfr(fit, "counts", ndraws = 50), "ggplot")
  expect_s3_class(pp_check_cfr(fit, "delay", ndraws = 50), "ggplot")

  # the observed death count should sit inside the posterior-predictive spread
  reps <- .cfr_replicate(fit, 200)
  pd <- reps$counts$n[reps$counts$outcome == "deaths"]
  obs <- reps$observed_counts$n[reps$observed_counts$outcome == "deaths"]
  expect_gte(obs, stats::quantile(pd, 0.01))
  expect_lte(obs, stats::quantile(pd, 0.99))
})

test_that("pp_check_cfr works when a formula covariate has missing values", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(30)
  ll <- simulate_linelist(n = 200, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  ll$age_group <- sample(c("young", "old"), nrow(ll), replace = TRUE)
  ll$age_group[c(5, 40)] <- NA
  d <- prepare_cfr_data(ll,
    obs_time = max(ll$onset_date) - 5, covariates = "age_group"
  )

  # brms drops the NA-covariate rows during fitting; pp_check_cfr() must still
  # be able to replay the real-time truncation on the rows it kept
  fit <- suppressWarnings(fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    formula = brms::bf(mu ~ 1, cfr ~ age_group),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, seed = 1
  ))

  expect_s3_class(pp_check_cfr(fit, "counts", ndraws = 50), "ggplot")
})
