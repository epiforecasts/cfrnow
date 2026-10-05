test_that("posterior_cumulative_incidence rejects objects that are not cfrnow fits", {
  expect_error(posterior_cumulative_incidence(1), "must come from fit_cfr")
})

test_that(".delay_cdf matches primarycensored at each time", {
  loc <- c(2.2, 2.6)
  sc <- c(0.5, 0.4)
  times <- c(0, 3, 10)
  f <- .delay_cdf("lognormal", times, 1, loc, sc)
  expect_equal(dim(f), c(2, 3))
  for (k in seq_along(times)) {
    expected <- primarycensored::pprimarycensored(
      rep(times[k], 2), stats::plnorm,
      pwindow = 1, L = -Inf, D = Inf, meanlog = loc, sdlog = sc,
      check = FALSE
    )
    expect_equal(f[, k], expected)
  }
})

test_that(".delay_cdf reaches 1 at a bounded delay's max", {
  f <- .delay_cdf("lognormal", c(10, 15, 20), 1, 2.4, 0.5, delay_max = 15)
  expected <- primarycensored::pprimarycensored(
    10, stats::plnorm,
    pwindow = 1, L = -Inf, D = 15, meanlog = 2.4, sdlog = 0.5
  )
  expect_equal(f[1, 1], expected)
  expect_equal(f[1, 2:3], c(1, 1))
})

test_that(".daily_hazards divides each day's incidence by the cases at risk", {
  eval_times <- 0:2
  inc <- list(
    death = array(c(0, 0.2, 0.3), c(1, 1, 3)),
    recovery = array(c(0, 0.1, 0.4), c(1, 1, 3))
  )
  cs <- .daily_hazards(inc, eval_times, 0:1, "hazard")
  expect_equal(as.vector(cs$death), c(0.2 / 1, 0.1 / 0.7))
  expect_equal(as.vector(cs$recovery), c(0.1 / 1, 0.3 / 0.7))

  sub <- .daily_hazards(inc, eval_times, 0:1, "subdistribution_hazard")
  expect_equal(as.vector(sub$death), c(0.2 / 1, 0.1 / 0.8))
  expect_equal(as.vector(sub$recovery), c(0.1 / 1, 0.3 / 0.9))
})

test_that(".daily_hazards gives NA once every case has resolved", {
  inc <- list(
    death = array(c(0.4, 0.4), c(1, 1, 2)),
    recovery = array(c(0.6, 0.6), c(1, 1, 2))
  )
  expect_true(is.na(.daily_hazards(inc, 0:1, 0, "hazard")$death[1]))
})

test_that(".covariate_patterns keeps one row per distinct pattern", {
  fit <- list(data = data.frame(
    y = 1:4, outcome = 0, pwindow = 1, swindow = 1,
    group = c("a", "b", "a", "b"), sex = c("f", "f", "f", "m")
  ))
  expect_equal(
    .covariate_patterns(fit),
    data.frame(group = c("a", "b", "b"), sex = c("f", "f", "m"))
  )
  no_covariates <- .covariate_patterns(list(data = fit$data[, 1:4]))
  expect_equal(dim(no_covariates), c(1, 0))
})

test_that(".with_addition_terms fills only the missing response variables", {
  nd <- .with_addition_terms(data.frame(group = "a", pwindow = 2))
  expect_equal(nd$pwindow, 2)
  expect_equal(nd$outcome, .CURE_CENSORED)
  expect_true(all(c("y", "swindow") %in% names(nd)))
})

test_that("posterior_cumulative_incidence is consistent with prob and the hazards", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(3)
  ll <- simulate_linelist(
    n = 400, cfr = 0.4, delay = LogNormal(2.4, 0.5),
    recovery = LogNormal(2.6, 0.4)
  )
  ll$group <- sample(c("a", "b"), nrow(ll), replace = TRUE)
  d <- prepare_cfr_data(ll, obs_time = NULL, covariates = "group")
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    recovery_delay = LogNormal(
      meanlog = Normal(2.6, 0.2), sdlog = Normal(0.4, 0.15)
    ),
    formula = brms::bf(mu ~ 1, prob ~ group),
    backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
    refresh = 0, seed = 1
  )
  groups <- data.frame(group = c("a", "b"))
  times <- 0:60

  ci <- posterior_cumulative_incidence(fit, groups, times = times)
  expect_setequal(names(ci), c(".draw", "row", "time", "outcome", "value"))
  ndraws <- posterior::ndraws(fit)
  expect_equal(nrow(ci), ndraws * 2 * length(times) * 2)

  death <- ci[ci$outcome == "death", ]
  recovery <- ci[ci$outcome == "recovery", ]
  # starts at 0, never decreases, and tends to prob (and 1 - prob) far out
  expect_true(all(death$value[death$time == 0] < 1e-6))
  ordered <- death[order(death$.draw, death$row, death$time), ]
  steps <- diff(ordered$value)[diff(ordered$time) == 1]
  expect_true(all(steps >= -1e-12))
  prob <- brms::posterior_linpred(fit,
    newdata = .with_addition_terms(groups), dpar = "prob", transform = TRUE
  )
  far_death <- death[death$time == 60, ]
  far_recovery <- recovery[recovery$time == 60, ]
  expect_equal(far_death$value, as.vector(prob), tolerance = 1e-3)
  expect_equal(far_death$value + far_recovery$value,
    rep(1, nrow(far_death)),
    tolerance = 1e-3
  )

  # the cause-specific hazards rebuild the cumulative incidence day by day
  h <- posterior_cumulative_incidence(fit, groups,
    times = 0:9, type = "hazard"
  )
  pick <- function(x, out, t) x$value[x$outcome == out & x$time == t]
  at_risk <- 1
  inc_death <- 0
  for (t in 0:9) {
    inc_death <- inc_death + at_risk * pick(h, "death", t)
    at_risk <- at_risk * (1 - pick(h, "death", t) - pick(h, "recovery", t))
    expect_equal(inc_death, pick(ci, "death", t + 1))
  }

  # by default, one curve per distinct covariate pattern in the fitted data
  expect_equal(
    sort(unique(posterior_cumulative_incidence(fit, times = 5)$row)), 1:2
  )

  # covariates only, no response columns, is enough for newdata
  expect_no_error(
    posterior_cumulative_incidence(fit, data.frame(group = "a"), times = 5)
  )
})

test_that("posterior_cumulative_incidence returns death only without a recovery delay", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(4)
  ll <- simulate_linelist(n = 200, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
    refresh = 0, seed = 1
  )
  nd <- fit$data[1, ]
  cs <- posterior_cumulative_incidence(fit, nd, times = 0:5, type = "hazard")
  sub <- posterior_cumulative_incidence(fit, nd,
    times = 0:5,
    type = "subdistribution_hazard"
  )
  expect_equal(unique(cs$outcome), "death")
  expect_equal(cs$value, sub$value)
})
