test_that("retrospective fit resolves every non-death and counts every death", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2, 3),
    death_date = as.Date(c("2026-01-10", NA, "2026-01-15", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = NULL)
  expect_equal(d$n_deaths, 2L)
  expect_equal(d$n_resolved, 2L) # both survivors treated as resolved
  expect_equal(d$n_cens, 0L) # nothing right-censored retrospectively
  expect_equal(d$n_cases, 4L)
})

test_that("real-time cut-off censors survivors and hides later deaths", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    # death 1 known by cut-off; death 2 happens after the cut-off; case 3 alive
    death_date = as.Date(c("2026-01-05", "2026-01-20", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  expect_equal(d$n_deaths, 1L) # only the death dated on/before cut-off
  expect_equal(d$n_cens, 2L) # the future death + the still-alive case
  expect_equal(d$n_resolved, 0L)
  expect_true(all(d$censor_time >= 0))
})

test_that("impossible onset->death delays are dropped with a warning", {
  ll <- data.frame(
    onset_date = as.Date(c("2026-01-10", "2026-01-01")),
    death_date = as.Date(c("2026-01-05", "2026-01-08")) # first: death before onset
  )
  expect_warning(d <- prepare_cfr_data(ll, obs_time = NULL), "unusable")
  expect_equal(d$n_deaths, 1L)
  expect_equal(d$n_cases, 1L) # the bad record is dropped entirely
})

test_that("missing onset dates are dropped with a warning, not silently kept", {
  ll <- data.frame(
    onset_date = as.Date(c("2026-01-01", NA, "2026-01-03")),
    death_date = as.Date(c("2026-01-10", "2026-01-08", NA))
  )
  expect_warning(d <- prepare_cfr_data(ll, obs_time = NULL), "unusable")
  expect_equal(d$n_cases, 2L) # the NA-onset row is dropped
  expect_equal(d$n_deaths, 1L) # not counted as a phantom death
  expect_false(anyNA(d$death_delay))
  expect_false(is.na(d$n_cases))
  expect_false(is.na(d$n_resolved))
})

test_that("inverted onset windows are dropped rather than producing negative widths", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-05"),
    onset_lower = as.Date("2026-01-05"),
    onset_upper = as.Date("2026-01-01"), # upper before lower
    death_date = as.Date("2026-01-10")
  )
  expect_warning(d <- prepare_cfr_data(ll, obs_time = NULL), "unusable")
  expect_equal(d$n_cases, 0L)
  expect_true(all(d$death_width >= 1))
})

test_that("cases with onset after the cut-off are excluded with a message", {
  ll <- data.frame(
    onset_date = as.Date(c("2026-01-01", "2026-01-20")), # second onsets post-cutoff
    death_date = as.Date(c("2026-01-08", NA))
  )
  expect_message(
    d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10")),
    "onset after the cut-off"
  )
  expect_equal(d$n_cases, 1L)
})

test_that("retrospective obs_time field is a Date, matching real-time mode", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01"),
    death_date = as.Date("2026-01-10")
  )
  expect_s3_class(prepare_cfr_data(ll, obs_time = NULL)$obs_time, "Date")
})

test_that("a non-NULL NA obs_time is rejected", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01"),
    death_date = as.Date("2026-01-10")
  )
  expect_error(prepare_cfr_data(ll, obs_time = NA), "valid date")
})

test_that("onset windows widen the primary censoring width", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01"),
    onset_lower = as.Date("2026-01-01"),
    onset_upper = as.Date("2026-01-04"),
    death_date = as.Date("2026-01-12")
  )
  d <- prepare_cfr_data(ll, obs_time = NULL)
  expect_equal(d$death_width, 4) # (upper - lower) + 1
})

test_that("recovered-by-cutoff cases are timed recoveries, not censored", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", NA, NA)),
    recovery_date = as.Date(c(NA, "2026-01-08", NA)) # case 2 recovered by cut-off
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  expect_equal(d$n_deaths, 1L) # case 1 died
  expect_equal(d$n_recovery, 1L) # case 2 recovered -> timed recovery
  expect_equal(d$recovery_delay, 6L) # onset 2026-01-02 to recovery 2026-01-08
  expect_equal(d$n_cens, 1L) # case 3 still unresolved -> censored
})

test_that("recovery_date absent reproduces the censor-everything behaviour", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1),
    death_date = as.Date(c(NA, NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  expect_equal(d$n_resolved, 0L) # no recovery info -> both censored
  expect_equal(d$n_cens, 2L)
})

test_that("a recovery before onset is dropped as unusable", {
  ll <- data.frame(
    onset_date = as.Date(c("2026-01-10", "2026-01-01")),
    death_date = as.Date(c(NA, NA)),
    recovery_date = as.Date(c("2026-01-05", "2026-01-09")) # case 1: recovers pre-onset
  )
  expect_warning(
    d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-15")), "unusable"
  )
  expect_equal(d$n_cases, 1L) # bad record dropped
  expect_equal(d$n_recovery, 1L) # case 2 is a valid recovery
})

test_that("onset and covariates are carried through to the cases frame", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", NA, NA)),
    region = c("A", "B", "A")
  )
  d <- prepare_cfr_data(ll,
    obs_time = as.Date("2026-01-10"),
    covariates = "region"
  )
  expect_true(all(c("y", "outcome", "pwindow", "swindow", "onset", "region")
  %in% names(d$cases)))
  expect_equal(nrow(d$cases), d$n_cases)
  expect_equal(sum(d$cases$outcome == 1), d$n_deaths)
  cure <- as_epidist_cure_model(d)
  expect_true("region" %in% names(cure))
  expect_s3_class(cure$onset, "Date")
})

test_that("a bounded delay simulates the recorded delays the model fits", {
  set.seed(9)
  mx <- 14
  ll <- simulate_linelist(
    n = 20000, cfr = 1, delay = LogNormal(2.4, 0.5, max = mx)
  )
  recorded <- as.numeric(ll$death_date - ll$onset_date)
  # the bound applies to the recorded delay, so it runs from 0 to max - 1
  expect_equal(max(recorded), mx - 1)

  # and each day carries the probability the likelihood gives it
  empirical <- as.numeric(table(factor(recorded, levels = 0:(mx - 1))))
  empirical <- empirical / sum(empirical)
  model <- primarycensored::dprimarycensored(
    0:(mx - 1), stats::plnorm,
    pwindow = 1, swindow = 1, D = mx, meanlog = 2.4, sdlog = 0.5
  )
  expect_lt(sum(abs(empirical - model)) / 2, 0.02)
})
