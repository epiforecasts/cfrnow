# Translate distspec delay / prob priors into the brms family and priors the
# epidist engine reads. distspec is cfrnow's user-facing idiom; this keeps the
# brms parameterisation an implementation detail.

`%||%` <- function(a, b) if (is.null(a)) b else a # nolint: coalesce_linter.

# A native delay parameter is either a fixed number or a Normal() prior.
.param_spec <- function(p) {
  if (is.numeric(p)) {
    return(list(fixed = p))
  }
  if (inherits(p, "dist_spec")) {
    dname <- get_distribution(p)
    pr <- get_parameters(p)
    if (dname == "fixed") {
      return(list(fixed = pr$value))
    }
    if (dname == "normal") {
      return(list(mean = pr$mean, sd = pr$sd))
    }
  }
  stop("delay parameters must be fixed numbers or Normal() priors.",
    call. = FALSE
  )
}

# One brms prior row for a dpar. `log_link` transforms to the link scale: a
# fixed value becomes a constant() prior, a Normal() becomes a Normal on the
# link scale (delta method for the sd when log-linked).
.prior_row <- function(spec, dpar, log_link) {
  tf <- function(x) if (log_link) log(x) else x
  prior_args <- if (nzchar(dpar)) {
    list(class = "Intercept", dpar = dpar)
  } else {
    list(class = "Intercept")
  }
  code <- if (!is.null(spec$fixed)) {
    sprintf("constant(%.8f)", tf(spec$fixed))
  } else {
    sc <- if (log_link) spec$sd / spec$mean else spec$sd
    sprintf("normal(%.6f, %.6f)", tf(spec$mean), sc)
  }
  do.call(brms::set_prior, c(list(code), prior_args))
}

# distspec Gamma() uses (shape, rate); brms Gamma uses (mu = mean, shape), both
# log-linked. Build the mu prior from mean = shape / rate.
.gamma_mu_prior <- function(sh, rate_spec, dpar) {
  if (!is.null(sh$fixed) && !is.null(rate_spec$fixed)) {
    return(.prior_row(list(fixed = sh$fixed / rate_spec$fixed), dpar, TRUE))
  }
  shc <- sh$fixed %||% sh$mean
  rtc <- rate_spec$fixed %||% rate_spec$mean
  vlog <- sqrt(((sh$sd %||% 0) / shc)^2 + ((rate_spec$sd %||% 0) / rtc)^2)
  prior_args <- if (nzchar(dpar)) {
    list(class = "Intercept", dpar = dpar)
  } else {
    list(class = "Intercept")
  }
  do.call(
    brms::set_prior,
    c(
      list(sprintf("normal(%.6f, %.6f)", log(shc / rtc), vlog)),
      prior_args
    )
  )
}

# distspec Weibull() uses (shape, scale); brms weibull uses (mu = mean, shape),
# both log-linked. Build the mu prior from mean = scale * gamma(1 + 1 / shape),
# propagating shape/scale uncertainty to log(mu) by the delta method.
.weibull_mu_prior <- function(shape_spec, scale_spec, dpar) {
  mu_of <- function(k, lambda) lambda * gamma(1 + 1 / k)
  if (!is.null(shape_spec$fixed) && !is.null(scale_spec$fixed)) {
    return(.prior_row(
      list(fixed = mu_of(shape_spec$fixed, scale_spec$fixed)), dpar, TRUE
    ))
  }
  kc <- shape_spec$fixed %||% shape_spec$mean
  lc <- scale_spec$fixed %||% scale_spec$mean
  # d log(mu)/d shape uses digamma(1 + 1/k); sign drops on squaring
  d_shape <- digamma(1 + 1 / kc) / kc^2
  vlog <- sqrt(
    ((scale_spec$sd %||% 0) / lc)^2 + (d_shape * (shape_spec$sd %||% 0))^2
  )
  prior_args <- if (nzchar(dpar)) {
    list(class = "Intercept", dpar = dpar)
  } else {
    list(class = "Intercept")
  }
  do.call(
    brms::set_prior,
    c(
      list(sprintf("normal(%.6f, %.6f)", log(lc) + lgamma(1 + 1 / kc), vlog)),
      prior_args
    )
  )
}

# A distspec delay -> list(family, prior, max). `main = TRUE` targets the death
# delay (mu is the main dpar); `main = FALSE` the recovery delay (r-prefixed
# dpars). `max` is the delay's upper bound, Inf when it has none.
.delay_family_prior <- function(delay, main = TRUE) {
  if (!inherits(delay, "dist_spec")) {
    stop("`delay` must be a distspec distribution, ",
      "e.g. LogNormal() or Gamma().",
      call. = FALSE
    )
  }
  fam <- get_distribution(delay)
  pars <- get_parameters(delay)
  loc_dpar <- if (main) "" else "rmu"
  scale_dpar <- if (main) "sigma" else "rsigma"
  shape_dpar <- if (main) "shape" else "rshape"
  out <- switch(fam,
    lognormal = {
      prior <- c(
        .prior_row(.param_spec(pars$meanlog), loc_dpar, FALSE),
        .prior_row(.param_spec(pars$sdlog), scale_dpar, TRUE)
      )
      list(family = brms::lognormal(), prior = prior)
    },
    gamma = {
      sh <- .param_spec(pars$shape)
      prior <- c(
        .gamma_mu_prior(sh, .param_spec(pars$rate), loc_dpar),
        .prior_row(sh, shape_dpar, TRUE)
      )
      list(family = stats::Gamma(link = "log"), prior = prior)
    },
    weibull = {
      sh <- .param_spec(pars$shape)
      prior <- c(
        .weibull_mu_prior(sh, .param_spec(pars$scale), loc_dpar),
        .prior_row(sh, shape_dpar, TRUE)
      )
      list(family = brms::weibull(), prior = prior)
    }
  )
  if (is.null(out)) {
    stop("cfrnow supports LogNormal(), Gamma() and Weibull() delays only.",
      call. = FALSE
    )
  }
  out$max <- .delay_max(delay)
  out
}

# The upper bound of a distspec delay, as a positive whole number of days, or
# Inf when the delay is unbounded. The bound is the longest delay the model can
# produce, so the likelihood truncates the delay distribution there. distspec
# rounds a fractional bound up when it discretises, so do the same here and
# keep one delay object to one support.
.delay_max <- function(delay) {
  # read the bound off the object: max() on a delay with priors resolves the
  # same value but messages about the uncertain parameters on the way
  mx <- attr(delay, "max")
  if (is.null(mx) || length(mx) != 1 || is.na(mx) || is.infinite(mx)) {
    return(Inf)
  }
  if (mx <= 0) {
    stop("a delay's `max` must be a positive number of days.", call. = FALSE)
  }
  ceiling(mx)
}

# A Normal(m, s) on logit(prob) whose induced mean/sd on the [0, 1] scale match
# `mean`/`sd` (moment-matched).
.prob_logitnormal <- function(prob_mean, prob_sd, class = "Intercept",
                              dpar = "prob") {
  moments <- function(m, s) {
    x <- seq(m - 6 * s, m + 6 * s, length.out = 2001)
    w <- stats::dnorm(x, m, s)
    w <- w / sum(w)
    p <- stats::plogis(x)
    pm <- sum(w * p)
    c(pm, sqrt(sum(w * (p - pm)^2)))
  }
  obj <- function(par) {
    sum((moments(par[1], exp(par[2])) - c(prob_mean, prob_sd))^2)
  }
  init <- c(
    stats::qlogis(prob_mean),
    log(prob_sd / (prob_mean * (1 - prob_mean) + 1e-6) + 1e-3)
  )
  opt <- stats::optim(init, obj, method = "Nelder-Mead")
  brms::set_prior(sprintf("normal(%.4f, %.4f)", opt$par[1], exp(opt$par[2])),
    class = class, dpar = dpar
  )
}

# Does the prob sub-formula carry a population intercept? A plain formula with
# no prob component (the default `mu ~ 1`) leaves prob intercept-only;
# `prob ~ 0 + x` or `prob ~ x - 1` drops it, in which case the prior belongs on
# the coefficients (class "b") rather than a non-existent intercept.
.prob_has_intercept <- function(formula) {
  prob_form <- if (inherits(formula, "brmsformula")) {
    formula$pforms$prob
  } else {
    NULL
  }
  if (is.null(prob_form)) {
    return(TRUE)
  }
  isTRUE(as.logical(attr(stats::terms(prob_form), "intercept")))
}

# A distspec Beta() prior -> the matching logit-scale brms prior. `class` is
# "Intercept" when the prob formula keeps its intercept, else "b" to place the
# prior on the (logit-scale) coefficients of an intercept-free formula.
.prob_prior_to_brms <- function(prob_prior, class = "Intercept",
                                dpar = "prob") {
  ok <- inherits(prob_prior, "dist_spec") &&
    get_distribution(prob_prior) == "beta"
  if (!ok) {
    stop("`", dpar, "_prior` must be a distspec Beta() distribution.",
      call. = FALSE
    )
  }
  p <- get_parameters(prob_prior)
  a <- p$shape1
  b <- p$shape2
  .prob_logitnormal(
    a / (a + b), sqrt(a * b / ((a + b)^2 * (a + b + 1))), class, dpar
  )
}

# Does `formula` carry an old-style `cfr ~ ...` sub-formula that needs
# translating to `prob ~ ...`? The warning is raised by the caller (fit_cfr()),
# not here, so lifecycle attributes it to the right environment.
.formula_has_cfr <- function(formula) {
  inherits(formula, "brmsformula") && !is.null(formula$pforms$cfr)
}

# Renames a `cfr ~ ...` sub-formula to `prob ~ ...`.
.rename_cfr_formula <- function(formula) {
  cfr_form <- formula$pforms$cfr
  cfr_form[[2]] <- quote(prob)
  formula$pforms$cfr <- NULL
  formula$pforms$prob <- cfr_form
  formula
}

# A loss-to-follow-up specification: for the death and the recovery half, the
# Stan expression for the probability of being lost (a dpar name or a literal),
# the dpars to estimate, and their priors. `NULL` means every outcome is
# eventually recorded; a single Beta() means one probability for both; a list
# gives each half its own fixed number or Beta().
.loss_spec <- function(loss_prior) {
  none <- list(death = "0", recovery = "0")
  if (is.null(loss_prior)) {
    return(list(
      parts = none, dpars = character(), prior = NULL, shared = FALSE
    ))
  }
  if (inherits(loss_prior, "dist_spec")) {
    # one probability whatever the outcome: the data identifies it
    return(list(
      parts = list(death = "loss", recovery = "loss"), dpars = "loss",
      prior = .prob_prior_to_brms(loss_prior, "Intercept", "loss"),
      shared = TRUE
    ))
  }
  if (!is.list(loss_prior) || !all(names(loss_prior) %in% names(none)) ||
        !all(names(none) %in% names(loss_prior))) {
    stop("`loss_prior` must be a Beta(), or a list with `death` and ",
      "`recovery` entries, each a Beta() or a number.",
      call. = FALSE
    )
  }
  dpar_of <- c(death = "dloss", recovery = "rloss")
  parts <- none
  dpars <- character()
  prior <- brms::empty_prior()
  for (nm in names(none)) {
    p <- loss_prior[[nm]]
    if (is.numeric(p)) {
      if (length(p) != 1 || is.na(p) || p < 0 || p >= 1) {
        stop("a fixed `loss_prior` entry must be a probability below 1.",
          call. = FALSE
        )
      }
      parts[[nm]] <- sprintf("%.8f", p)
    } else {
      parts[[nm]] <- dpar_of[[nm]]
      dpars <- c(dpars, dpar_of[[nm]])
      prior <- c(
        prior, .prob_prior_to_brms(p, "Intercept", dpar_of[[nm]])
      )
    }
  }
  list(parts = parts, dpars = dpars, prior = prior, shared = FALSE)
}

# Is this half of the loss specification switched off (a fixed zero)?
.loss_is_off <- function(part) {
  !is.na(suppressWarnings(as.numeric(part))) && as.numeric(part) == 0
}
