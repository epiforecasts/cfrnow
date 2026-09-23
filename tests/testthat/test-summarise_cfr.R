test_that(".ascertainment_adjust shifts the CFR the documented way", {
  cfr <- c(0.05, 0.1, 0.4)
  logit <- stats::qlogis(cfr)

  # r = 1 is a no-op
  expect_equal(.ascertainment_adjust(logit, 1), cfr)

  # r > 1 (fatal cases over-ascertained) lowers the CFR; r < 1 raises it
  expect_true(all(.ascertainment_adjust(logit, 2) < cfr))
  expect_true(all(.ascertainment_adjust(logit, 0.5) > cfr))

  # it is exactly a -log(r) shift on the logit scale
  expect_equal(
    .ascertainment_adjust(logit, 2),
    stats::plogis(logit - log(2))
  )

  # inverse ratios are symmetric on the logit scale
  expect_equal(
    stats::qlogis(.ascertainment_adjust(logit, 3)) +
      stats::qlogis(.ascertainment_adjust(logit, 1 / 3)),
    2 * logit
  )
})

test_that("summary() rejects an invalid ascertainment_ratio", {
  stub <- structure(list(cfrnow = list()), class = "cfrnow_fit")
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(
      summary(stub, ascertainment_ratio = bad),
      "positive number"
    )
  }
})

test_that("summary() applies the ascertainment correction to a real fit", {
  testthat::skip_if_not_installed("cmdstanr")
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

  base <- summary(fit)
  down <- summary(fit, ascertainment_ratio = 2)
  up <- summary(fit, ascertainment_ratio = 0.5)

  cfr_med <- function(s) s[s$quantity == "prob", "q50"]
  expect_lt(cfr_med(down), cfr_med(base))
  expect_gt(cfr_med(up), cfr_med(base))
  # the default leaves the CFR untouched
  expect_equal(cfr_med(summary(fit, ascertainment_ratio = 1)), cfr_med(base))
  expect_equal(attr(down, "ascertainment_ratio"), 2)

  # the weak-identification flag reflects the fit, not the ascertainment lens
  expect_equal(
    attr(down, "cfr_low_information"),
    attr(base, "cfr_low_information")
  )
})

test_that(".prob_is_grouped detects a grouped prob formula", {
  mk <- function(f) structure(list(formula = f), class = "brmsfit")
  # prob defaults to intercept-only
  expect_false(.prob_is_grouped(mk(brms::bf(mu ~ 1))))
  expect_false(.prob_is_grouped(mk(brms::bf(mu ~ 1, prob ~ 1))))
  expect_true(.prob_is_grouped(mk(brms::bf(mu ~ 1, prob ~ site))))
  expect_true(.prob_is_grouped(mk(brms::bf(mu ~ 1, prob ~ 0 + site))))
  expect_true(.prob_is_grouped(mk(brms::bf(mu ~ 1, prob ~ (1 | site)))))
})
