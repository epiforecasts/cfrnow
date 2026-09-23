# Input-validation guards and pure-R helpers. None of these need CmdStan, so
# they run everywhere and exercise the error paths taken before the sampler.

test_that("as_epidist_cure_model needs the model columns", {
  expect_error(
    as_epidist_cure_model(data.frame(x = 1)),
    "y|outcome|pwindow|swindow"
  )
})

test_that("as_epidist_cure_model builds cure rows from prepare_cfr_data output", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", "2026-01-20", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  cure <- as_epidist_cure_model(d)
  expect_s3_class(cure, "epidist_cure_model")
  expect_equal(sum(cure$outcome == 1), d$n_deaths) # 1 observed death
  expect_equal(sum(cure$outcome == 0), d$n_cens) # future death + survivor
  expect_equal(nrow(cure), d$n_cases)
})

test_that("recovery dates produce timed-recovery rows flagged for two outcomes", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", NA, NA)),
    recovery_date = as.Date(c(NA, "2026-01-08", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  cure <- as_epidist_cure_model(d)
  expect_true(attr(cure, "use_recovery"))
  expect_equal(sum(cure$outcome == 2), d$n_recovery) # timed recovery row
  expect_equal(sum(cure$outcome == 1), d$n_deaths)
})

test_that("naive_cfr is deaths/cases, and NA when there are no cases", {
  expect_equal(naive_cfr(3, 10), 0.3)
  expect_true(is.na(naive_cfr(0, 0)))
})

test_that(".prob_prior_sd reads a normal prob prior and is NA otherwise", {
  expect_true(is.na(.prob_prior_sd(brms::set_prior("beta(1, 1)", class = "b"))))
  s <- .prob_prior_sd(.prob_prior_to_brms(Beta(1, 1)))
  expect_true(is.finite(s) && s > 0 && s < 0.5) # (0,1)-scale prior sd
})

test_that(".delay_family_prior maps fixed params to constant priors", {
  p <- .delay_family_prior(LogNormal(meanlog = 2.41, sdlog = 0.51))$prior
  expect_true(all(grepl("^constant\\(", p$prior)))
  loc <- as.numeric(sub("constant\\(([^)]+)\\)", "\\1", p$prior[p$dpar == ""]))
  expect_equal(loc, 2.41, tolerance = 1e-6) # meanlog is identity-linked
})

test_that(".delay_family_prior maps Normal() params to normal priors and family", {
  dd <- .delay_family_prior(LogNormal(
    meanlog = Normal(2.4, 0.2),
    sdlog = Normal(0.5, 0.15)
  ))
  expect_equal(dd$family$family, "lognormal")
  expect_true(any(grepl("normal\\(2.4", dd$prior$prior)))
  expect_error(.delay_family_prior(5), "distspec")
})

test_that(".prob_prior_to_brms turns a Beta into a logit-normal of matching sd", {
  p <- .prob_prior_to_brms(Beta(1, 1)) # uniform: sd = 1/sqrt(12)
  expect_equal(.prob_prior_sd(p), 1 / sqrt(12), tolerance = 0.03)
  expect_error(.prob_prior_to_brms(Gamma(1, 1)), "Beta")
})

test_that("simulate_linelist requires a delay", {
  expect_error(simulate_linelist(n = 5), "supply a `delay`")
})

test_that("prepare_cfr_data requires onset_date and death_date columns", {
  expect_error(
    prepare_cfr_data(data.frame(death_date = as.Date("2026-01-05"))),
    "onset_date"
  )
  expect_error(
    prepare_cfr_data(data.frame(onset_date = as.Date("2026-01-01"))),
    "death_date"
  )
})

test_that("summary.cfrnow_fit rejects objects not from fit_cfr", {
  expect_error(summary.cfrnow_fit(1), "fit_cfr")
})

test_that(".assert_within_max rejects data outside a bounded delay", {
  cure <- as_epidist_cure_model(data.frame(
    y = c(5L, 40L, 10L, 90L),
    outcome = c(.CURE_DEATH, .CURE_DEATH, .CURE_RECOVERY, .CURE_CENSORED),
    pwindow = 1, swindow = 1
  ))
  expect_true(.assert_within_max(cure, Inf, Inf))
  expect_error(.assert_within_max(cure, 30, Inf), "1 death\\(s\\)")
  expect_error(.assert_within_max(cure, Inf, 5), "1 recovery\\(ies\\)")

  # a delay whose secondary window closes past the bound is rejected too: the
  # 40-day death needs a max of 41, and a wider secondary window needs more
  expect_true(.assert_within_max(cure, 41, 95))
  expect_error(.assert_within_max(cure, 40, 95), "1 death\\(s\\)")
  wide <- cure
  wide$swindow <- 2
  expect_error(.assert_within_max(wide, 41, 95), "1 death\\(s\\)")

  # a censored case beyond every bound cannot resolve in a two-outcome fit
  expect_error(.assert_within_max(cure, 60, 60), "still unresolved")
  # death-only: an unresolved case is simply a survivor, whatever its follow-up
  attr(cure, "use_recovery") <- FALSE
  expect_true(.assert_within_max(cure, 60, NULL))
})

test_that(".assert_loss_identified needs unresolved cases", {
  resolved <- as_epidist_cure_model(data.frame(
    y = c(5L, 0L), outcome = c(.CURE_DEATH, .CURE_RESOLVED),
    pwindow = 1, swindow = 1
  ))
  expect_error(
    .assert_loss_identified(resolved, use_recovery = TRUE), "still unresolved"
  )

  censored <- as_epidist_cure_model(data.frame(
    y = c(5L, 20L), outcome = c(.CURE_DEATH, .CURE_CENSORED),
    pwindow = 1, swindow = 1
  ))
  expect_true(.assert_loss_identified(censored, use_recovery = TRUE))
  expect_warning(
    .assert_loss_identified(censored, use_recovery = FALSE),
    "weakly identified"
  )
})
