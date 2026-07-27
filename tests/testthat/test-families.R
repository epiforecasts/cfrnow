test_that("simulated delays match the requested mean and sd", {
  set.seed(1)
  specs <- list(
    LogNormal(mean = 12.75, sd = 7),
    Gamma(mean = 12.75, sd = 7),
    Weibull(mean = 12.75, sd = 7)
  )
  for (delay in specs) {
    ll <- simulate_linelist(n = 20000, cfr = 1, delay = delay)
    d <- as.numeric(ll$death_date - ll$onset_date)
    expect_equal(mean(d), 12.75, tolerance = 0.3)
    expect_equal(sd(d), 7, tolerance = 0.3)
  }
})

test_that("a recovery delay adds recoveries for non-fatal cases only", {
  set.seed(1)
  ll <- simulate_linelist(
    n = 300, cfr = 0.5,
    delay = LogNormal(mean = 12.75, sd = 7),
    recovery = LogNormal(mean = 21, sd = 9)
  )
  expect_true("recovery_date" %in% names(ll))
  fatal <- !is.na(ll$death_date)
  expect_true(all(is.na(ll$recovery_date[fatal]))) # fatal: no recovery
  expect_true(all(!is.na(ll$recovery_date[!fatal]))) # non-fatal: recovered
  no_rec <- simulate_linelist(n = 5, delay = LogNormal(2, 0.5))
  expect_false("recovery_date" %in% names(no_rec))
})

test_that("unsupported families and prior-parameter delays are rejected", {
  expect_error(sample_delay(5, Exponential(rate = 1)), "supports") # exp not supported
  expect_error(sample_delay(5, LogNormal(
    meanlog = Normal(2, 0.1), sdlog = Normal(0.5, 0.1)
  )), "fixed")
})
