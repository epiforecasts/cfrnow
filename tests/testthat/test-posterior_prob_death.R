test_that("posterior_prob_death rejects objects that are not cfrnow fits", {
  expect_error(posterior_prob_death(1), "must come from fit_cfr")
  expect_error(
    posterior_prob_death(structure(list(), class = "brmsfit")),
    "must come from fit_cfr"
  )
})

test_that(".delay_survivor matches primarycensored for each supported family", {
  nd <- 4
  y <- c(3, 8)
  pwindow <- c(1, 2)

  loc <- matrix(c(2.2, 2.4, 2.3, 2.5, 2.1, 2.6, 2.0, 2.7), nd, 2)
  sc <- matrix(0.5, nd, 2)
  s <- .delay_survivor("lognormal", y, pwindow, loc, sc)
  for (j in 1:2) {
    expected <- 1 - primarycensored::pprimarycensored(
      rep(y[j], nd), stats::plnorm,
      pwindow = pwindow[j], L = -Inf, D = Inf,
      meanlog = loc[, j], sdlog = sc[, j], check = FALSE
    )
    expect_equal(s[, j], expected)
  }

  loc_g <- matrix(c(8, 9, 10, 11, 7, 12, 6, 13), nd, 2)
  shape_g <- matrix(3, nd, 2)
  s_g <- .delay_survivor("gamma", y, pwindow, loc_g, shape_g)
  for (j in 1:2) {
    expected <- 1 - primarycensored::pprimarycensored(
      rep(y[j], nd), stats::pgamma,
      pwindow = pwindow[j], L = -Inf, D = Inf,
      shape = shape_g[, j], rate = shape_g[, j] / loc_g[, j], check = FALSE
    )
    expect_equal(s_g[, j], expected)
  }

  loc_w <- matrix(c(8, 9, 10, 11, 7, 12, 6, 13), nd, 2)
  shape_w <- matrix(1.5, nd, 2)
  s_w <- .delay_survivor("weibull", y, pwindow, loc_w, shape_w)
  for (j in 1:2) {
    scale_w <- loc_w[, j] / gamma(1 + 1 / shape_w[, j])
    expected <- 1 - primarycensored::pprimarycensored(
      rep(y[j], nd), stats::pweibull,
      pwindow = pwindow[j], L = -Inf, D = Inf,
      shape = shape_w[, j], scale = scale_w, check = FALSE
    )
    expect_equal(s_w[, j], expected)
  }
})

test_that(".delay_survivor decreases in y (a longer follow-up survives less)", {
  loc <- matrix(2.4, 5, 1)
  sc <- matrix(0.5, 5, 1)
  s_short <- .delay_survivor("lognormal", 2, 1, loc, sc)
  s_long <- .delay_survivor("lognormal", 20, 1, loc, sc)
  expect_true(all(s_long < s_short))
})

test_that(".pi_death applies Bayes' rule to censored cases only", {
  nd <- 3
  n <- 4
  prob <- matrix(0.3, nd, n)
  s_d <- matrix(0.8, nd, n)
  s_r <- matrix(0.6, nd, n)
  outcome <- c(.CURE_DEATH, .CURE_RECOVERY, .CURE_RESOLVED, .CURE_CENSORED)

  pi <- .pi_death(prob, s_d, s_r, outcome)
  expect_true(all(pi[, 1] == 1)) # observed death
  expect_true(all(pi[, 2] == 0)) # observed recovery
  expect_true(all(pi[, 3] == 0)) # resolved non-death

  expected_censored <- 0.3 * 0.8 / (0.3 * 0.8 + 0.7 * 0.6)
  expect_equal(pi[, 4], rep(expected_censored, nd))
})

test_that(".pi_death reduces to S_R = 1 when there is no recovery delay", {
  prob <- matrix(0.4, 2, 1)
  s_d <- matrix(0.9, 2, 1)
  pi <- .pi_death(prob, s_d, 1, .CURE_CENSORED)
  expect_equal(pi[, 1], rep(0.4 * 0.9 / (0.4 * 0.9 + 0.6), 2))
})

test_that("posterior_prob_death gives 0/1 draws matching outcome on a fully resolved fit", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(1)
  ll <- simulate_linelist(n = 300, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL) # every case resolved: no censoring
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, seed = 1
  )

  pi <- posterior_prob_death(fit)
  pi_df <- as.data.frame(pi)
  pi_mat <- unname(as.matrix(pi_df[, grep("^pi\\[", names(pi_df))]))
  expected <- as.numeric(fit$data$outcome == .CURE_DEATH)
  # every draw reproduces the observed outcome exactly (no censored cases)
  expect_true(all(apply(pi_mat, 1, identical, expected)))
  # averaging over cases reproduces the naive CFR for every draw
  expect_equal(rowMeans(pi_mat), rep(mean(expected), nrow(pi_mat)))
})

test_that("posterior_prob_death gives an unresolved case's pi near prob or 0/1", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(2)
  ll <- simulate_linelist(
    n = 600, cfr = 0.4, onset_days = 40, delay = LogNormal(2.4, 0.5)
  )
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, seed = 1
  )

  pi <- posterior_prob_death(fit)
  pi_df <- as.data.frame(pi)
  pi_mat <- as.matrix(pi_df[, grep("^pi\\[", names(pi_df))])
  cens <- which(fit$data$outcome == .CURE_CENSORED)
  expect_gt(length(cens), 0)

  # a censored case's pi lies strictly inside (0, 1)
  expect_true(all(pi_mat[, cens] > 0 & pi_mat[, cens] < 1))

  # the least-followed-up censored case (smallest y) sits closest to prob;
  # its pi should not have moved as far towards 0 as the longest-followed one
  y_cens <- fit$data$y[cens]
  newest <- cens[which.min(y_cens)]
  oldest <- cens[which.max(y_cens)]
  prob_draws <- brms::posterior_linpred(
    fit,
    dpar = "prob", transform = TRUE
  )
  dist_newest <- abs(rowMeans(prob_draws) - rowMeans(pi_mat[, newest, drop = FALSE]))
  dist_oldest <- abs(rowMeans(prob_draws) - rowMeans(pi_mat[, oldest, drop = FALSE]))
  expect_lt(mean(dist_newest), mean(dist_oldest))
})
