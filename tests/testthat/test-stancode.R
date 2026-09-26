# Coverage without a sampler: drive model generation through make_stancode,
# which exercises the family / formula / stancode / translate paths (cure_model.R
# and translate.R) but needs no CmdStan, so it runs on CI where the fit tests skip.

stancode_for <- function(cure, delay, recovery_delay = NULL, extra = NULL) {
  dd <- .delay_family_prior(delay)
  prior <- c(.prob_prior_to_brms(Beta(1, 1)), dd$prior, extra)
  if (!is.null(recovery_delay) && isTRUE(attr(cure, "use_recovery"))) {
    rd <- .delay_family_prior(recovery_delay, main = FALSE)
    prior <- c(prior, rd$prior)
    attr(cure, "recovery_family") <- brms:::validate_family(rd$family)
  }
  epidist::epidist(cure,
    formula = mu ~ 1, family = dd$family, prior = prior,
    merge_priors = FALSE, fn = brms::make_stancode
  )
}

test_that("a death-only lognormal fit generates the cure lpmf", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(n = 100, cfr = 0.4, delay = LogNormal(2.4, 0.5)),
    obs_time = NULL
  ))
  code <- stancode_for(cure, LogNormal(meanlog = 2.41, sdlog = 0.51))
  expect_true(grepl("cfrnow_lognormal_lpmf", code))
  expect_true(grepl("primarycensored", code))
})

test_that("a gamma fit generates a gamma cure lpmf", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(n = 100, cfr = 0.4, delay = Gamma(mean = 8, sd = 4)),
    obs_time = NULL
  ))
  code <- stancode_for(cure, Gamma(shape = Normal(4, 1), rate = Normal(0.5, 0.2)))
  expect_true(grepl("cfrnow_gamma_lpmf", code))
})

test_that("a two-outcome fit generates recovery branches with its own family", {
  ll <- simulate_linelist(
    n = 200, cfr = 0.4, delay = Gamma(mean = 12, sd = 6),
    recovery = LogNormal(mean = 20, sd = 8)
  )
  cure <- as_epidist_cure_model(prepare_cfr_data(ll, obs_time = NULL))
  expect_true(attr(cure, "use_recovery"))
  code <- stancode_for(
    cure, Gamma(shape = Normal(4, 1), rate = Normal(0.3, 0.1)),
    recovery_delay = LogNormal(meanlog = Normal(2.9, 0.3), sdlog = Normal(0.5, 0.2))
  )
  expect_true(grepl("outcome == 2", code)) # timed-recovery branch present
  expect_true(grepl("rmu", code)) # recovery params are r-prefixed
})

test_that("an intercept-free prob formula routes the prior to the coefficients", {
  # where prob_prior lands depends on whether the prob formula keeps its
  # intercept
  expect_true(.prob_has_intercept(mu ~ 1)) # prob defaults to intercept-only
  expect_true(.prob_has_intercept(brms::bf(mu ~ 1, prob ~ grp)))
  expect_false(.prob_has_intercept(brms::bf(mu ~ 1, prob ~ 0 + grp)))

  expect_equal(.prob_prior_to_brms(Beta(1, 1))$class, "Intercept")
  expect_equal(.prob_prior_to_brms(Beta(1, 1), "b")$class, "b")

  a <- simulate_linelist(n = 80, cfr = 0.3, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 80, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  ca <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  cb <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  ca$grp <- "x"
  cb$grp <- "y"
  cure <- as_epidist_cure_model(rbind(ca, cb))

  dd <- .delay_family_prior(LogNormal(meanlog = 2.41, sdlog = 0.51))
  f <- brms::bf(mu ~ 1, prob ~ 0 + grp)

  # the prior on the (absent) intercept is what brms rejects ...
  expect_error(
    epidist::epidist(cure,
      formula = f, family = dd$family,
      prior = c(.prob_prior_to_brms(Beta(1, 1), "Intercept"), dd$prior),
      merge_priors = FALSE, fn = brms::make_stancode
    )
  )
  # ... and moving it onto the coefficients generates cleanly
  code <- epidist::epidist(cure,
    formula = f, family = dd$family,
    prior = c(.prob_prior_to_brms(Beta(1, 1), "b"), dd$prior),
    merge_priors = FALSE, fn = brms::make_stancode
  )
  expect_true(grepl("prob", code))
})

test_that(".formula_has_cfr / .rename_cfr_formula translate `cfr ~ ...`", {
  expect_false(.formula_has_cfr(mu ~ 1))
  expect_false(.formula_has_cfr(brms::bf(mu ~ 1, prob ~ grp)))
  expect_true(.formula_has_cfr(brms::bf(mu ~ 1, cfr ~ grp)))

  translated <- .rename_cfr_formula(brms::bf(mu ~ 1, cfr ~ 0 + grp))
  expect_null(translated$pforms$cfr)
  expect_equal(translated$pforms$prob, prob ~ 0 + grp)
})

test_that("a delay max truncates the generated lpmf", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(n = 100, cfr = 0.4, delay = LogNormal(2.4, 0.5)),
    obs_time = NULL
  ))
  # primarycensored's own functions mention positive_infinity(), so look only
  # at the cure lpmf cfrnow generates
  cure_lpmf <- function(code) {
    regmatches(code, regexpr(
      "(?s)real cfrnow_lognormal_lpmf.*?\n}", code,
      perl = TRUE
    ))
  }
  unbounded <- cure_lpmf(stancode_for(
    cure, LogNormal(meanlog = 2.41, sdlog = 0.51)
  ))
  expect_true(grepl("positive_infinity()", unbounded, fixed = TRUE))

  attr(cure, "delay_max") <- 30
  bounded <- cure_lpmf(stancode_for(
    cure, LogNormal(meanlog = 2.41, sdlog = 0.51, max = 30)
  ))
  # distspec truncates at max, so recorded delays run from 0 to max - 1
  expect_true(grepl("30.00000000", bounded, fixed = TRUE))
  expect_false(grepl("positive_infinity()", bounded, fixed = TRUE))
})

test_that(".delay_max reads a bound without resolving uncertain parameters", {
  expect_identical(.delay_max(LogNormal(2.4, 0.5)), Inf)
  expect_identical(.delay_max(LogNormal(2.4, 0.5, max = 30)), 30)
  expect_identical(
    .delay_max(Gamma(shape = Normal(2, 1), rate = Normal(0.3, 0.1), max = 14)),
    14
  )
  expect_error(.delay_max(LogNormal(2.4, 0.5, max = 0)), "positive number")

  # distspec rounds a fractional bound up when it discretises, so the fit has
  # to agree or the same object means two supports
  expect_identical(.delay_max(LogNormal(2.4, 0.5, max = 34.5)), 35)
  expect_identical(
    length(get_pmf(discretise(LogNormal(2.4, 0.5, max = 34.5)))), 35L
  )
})

test_that("survival terms are computed on the log scale", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(
      n = 200, cfr = 0.4, delay = Gamma(mean = 6, sd = 5),
      recovery = Gamma(mean = 12, sd = 4)
    ),
    obs_time = as.Date("2026-02-01")
  ))
  code <- stancode_for(
    cure, Gamma(shape = Normal(1.4, 0.5), rate = Normal(0.2, 0.1)),
    Gamma(shape = Normal(8, 3), rate = Normal(0.7, 0.3))
  )
  # a long follow-up rounds the cdf to 1, where log1m() rejects the draw
  expect_true(grepl("log1m_exp", code, fixed = TRUE))
  expect_true(grepl("primarycensored_lcdf", code, fixed = TRUE))
  expect_false(grepl("log1m(fbar", code, fixed = TRUE))
})

test_that("a loss_prior adds loss parameters and the mixture to the lpmf", {
  ll <- simulate_linelist(
    n = 100, cfr = 0.4, delay = Gamma(mean = 6, sd = 4),
    recovery = Gamma(mean = 12, sd = 4)
  )
  cure <- as_epidist_cure_model(
    prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
  )
  code_for <- function(loss_prior) {
    spec <- .loss_spec(loss_prior)
    attr(cure, "loss_parts") <- spec$parts
    attr(cure, "loss_dpars") <- spec$dpars
    stancode_for(
      cure, Gamma(shape = Normal(2, 1), rate = Normal(0.3, 0.1)),
      Gamma(shape = Normal(9, 3), rate = Normal(0.75, 0.3)), spec$prior
    )
  }

  plain <- code_for(NULL)
  expect_false(grepl("real loss", plain, fixed = TRUE))

  # one probability whatever the outcome
  shared <- code_for(Beta(1, 1))
  expect_true(grepl("real loss,", shared, fixed = TRUE))
  expect_true(grepl("log1m(loss) + log(prob)", shared, fixed = TRUE))
  expect_true(grepl("log(prob) + log(loss)", shared, fixed = TRUE))
  expect_true(grepl("log1m(prob) + log(loss)", shared, fixed = TRUE))

  # every death recorded, a recovery sometimes not: no loss term for deaths
  one_sided <- code_for(list(death = 0, recovery = Beta(1, 1)))
  expect_false(grepl("dloss", one_sided, fixed = TRUE))
  expect_true(grepl("log1m(rloss) + log1m(prob)", one_sided, fixed = TRUE))
  expect_true(grepl("log(prob) + log_surv_d", one_sided, fixed = TRUE))
  expect_false(grepl("log(prob) + log(", one_sided, fixed = TRUE))

  # a fixed probability becomes a literal rather than a parameter
  fixed <- code_for(list(death = 0.05, recovery = Beta(1, 1)))
  expect_false(grepl("real dloss", fixed, fixed = TRUE))
  expect_true(grepl("log1m(0.05000000)", fixed, fixed = TRUE))
})

test_that(".loss_spec rejects malformed loss priors", {
  expect_error(.loss_spec(list(death = 0)), "`death` and `recovery`")
  expect_error(.loss_spec(list(death = 0, recovery = 1)), "probability below 1")
  expect_error(
    .loss_spec(list(death = 0, recovery = -0.1)), "probability below 1"
  )
})
