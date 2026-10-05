#' Posterior cumulative incidence and hazards of death and recovery
#'
#' Turns a fit's outcome probability and delay distributions into the
#' quantities a competing-risks analysis reports, for each covariate pattern
#' in `newdata`. The cumulative incidence of death and of recovery by time `t`
#' since onset are
#' \deqn{I_D(t) = p \, F_D(t), \qquad I_R(t) = (1 - p) \, F_R(t)}
#' where `p` is `prob` and `F_D` and `F_R` are the distribution functions of
#' the death and recovery delays, censored for the day of onset and truncated
#' at each delay's `max` as in the likelihood.
#'
#' The hazards are daily, matching the daily resolution of the data. The
#' cause-specific hazard of death on day `t` is the probability of dying
#' during `[t, t + 1)` among cases still unresolved at `t`,
#' \deqn{h_D(t) = \frac{I_D(t + 1) - I_D(t)}{1 - I_D(t) - I_R(t)},}
#' and the subdistribution (Fine-Gray) hazard keeps cases that have recovered
#' in the risk set,
#' \deqn{\tilde h_D(t) = \frac{I_D(t + 1) - I_D(t)}{1 - I_D(t)}.}
#' Recovery works the same way. Hazard ratios between two covariate patterns
#' are ratios of these draws. They vary with `t`, since the model does not
#' assume proportional hazards; covariates act on `prob` and on the delays
#' separately.
#'
#' A fit without a `recovery_delay` returns death only. Its survivors never
#' resolve in the model, so the cause-specific and subdistribution hazards
#' coincide. Loss to follow-up only changes which outcomes get recorded and
#' plays no part here.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @param newdata A data frame with the covariates of the `prob` and delay
#'   formulas, one row per covariate pattern. Defaults to the fitted cases,
#'   `object$data`.
#' @param times Days since onset at which to evaluate. Defaults to every day
#'   from 0 to the longest follow-up in the fitted data.
#' @param type `"cumulative_incidence"` (the default), `"hazard"` for the
#'   cause-specific hazard or `"subdistribution_hazard"`.
#' @param pwindow Width of the onset window, in days, that the delays are
#'   censored over. Defaults to a single day.
#' @return A data frame with one row per posterior draw, row of `newdata`,
#'   time and outcome: `.draw`, `row` (the row of `newdata`), `time`,
#'   `outcome` (`"death"` or `"recovery"`) and `value`. A hazard is `NA` at a
#'   time when every case has already resolved.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' ll$group <- sample(c("a", "b"), nrow(ll), replace = TRUE)
#' d <- prepare_cfr_data(ll,
#'   obs_time = max(ll$onset_date) - 5,
#'   covariates = "group"
#' )
#' fit <- fit_cfr(d,
#'   delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)),
#'   formula = brms::bf(mu ~ group, prob ~ group)
#' )
#' groups <- data.frame(group = c("a", "b"))
#'
#' # cumulative incidence of death by day 7, 14 and 28 in each group
#' ci <- posterior_cumulative_incidence(fit, groups, times = c(7, 14, 28))
#' aggregate(value ~ row + time, ci, stats::median)
#'
#' # cause-specific hazard ratio of death, group b vs a, over time
#' h <- posterior_cumulative_incidence(fit, groups,
#'   times = 0:28,
#'   type = "hazard"
#' )
#' a <- h[h$outcome == "death" & h$row == 1, ]
#' b <- h[h$outcome == "death" & h$row == 2, ]
#' a$hr <- b$value / a$value
#' aggregate(hr ~ time, a, stats::median)
#' }
#' @family fit
#' @export
posterior_cumulative_incidence <- function(object, newdata = NULL,
                                           times = NULL,
                                           type = c(
                                             "cumulative_incidence",
                                             "hazard",
                                             "subdistribution_hazard"
                                           ),
                                           pwindow = 1) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  type <- match.arg(type)
  newdata <- .with_addition_terms(newdata %||% object$data)
  times <- times %||% seq(0, ceiling(max(object$data$y)))
  checkmate::assert_numeric(times, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_number(pwindow, lower = 0)

  lp <- function(dpar) {
    brms::posterior_linpred(object,
      newdata = newdata, dpar = dpar,
      transform = TRUE
    )
  }
  # a hazard on day t needs the cumulative incidence at t + 1 as well
  eval_times <- if (type == "cumulative_incidence") {
    times
  } else {
    sort(unique(c(times, times + 1)))
  }
  prob <- lp("prob")
  incidence <- list(death = .cumulative_incidence(
    prob, object$cfrnow$family, c("mu", "sigma", "shape"),
    lp, eval_times, pwindow, object$cfrnow$delay_max %||% Inf
  ))
  if (isTRUE(object$cfrnow$use_recovery)) {
    incidence$recovery <- .cumulative_incidence(
      1 - prob, object$cfrnow$recovery_family, c("rmu", "rsigma", "rshape"),
      lp, eval_times, pwindow, object$cfrnow$recovery_max %||% Inf
    )
  }

  values <- if (type == "cumulative_incidence") {
    now <- match(times, eval_times)
    lapply(incidence, function(inc) inc[, , now, drop = FALSE])
  } else {
    .daily_hazards(incidence, eval_times, times, type)
  }

  draws <- posterior::as_draws_df(object)$.draw
  out <- lapply(names(values), function(outcome) {
    v <- values[[outcome]]
    data.frame(
      .draw = rep(draws, times = dim(v)[2] * dim(v)[3]),
      row = rep(rep(seq_len(dim(v)[2]), each = dim(v)[1]), times = dim(v)[3]),
      time = rep(times, each = dim(v)[1] * dim(v)[2]),
      outcome = outcome,
      value = as.vector(v)
    )
  })
  do.call(rbind, out)
}

# brms checks newdata for the variables in the response's addition terms even
# though the linear predictors do not use them, so fill in placeholders for a
# newdata that holds covariates only.
.with_addition_terms <- function(newdata) {
  placeholders <- list(
    y = 0, outcome = .CURE_CENSORED, pwindow = 1, swindow = 1
  )
  for (v in setdiff(names(placeholders), names(newdata))) {
    newdata[[v]] <- placeholders[[v]]
  }
  newdata
}

# Cumulative incidence of one outcome: its probability (an ndraws x nrows
# matrix) times its delay's distribution function at each time in `eval_times`.
# `dpars` names the delay's location and its lognormal / other scale
# parameter. Returns an ndraws x nrows x length(eval_times) array.
.cumulative_incidence <- function(prob, family, dpars, lp, eval_times, pwindow,
                                  delay_max) {
  loc <- lp(dpars[1])
  shape_or_sd <- lp(if (family == "lognormal") dpars[2] else dpars[3])
  inc <- array(NA_real_, c(nrow(prob), ncol(prob), length(eval_times)))
  for (j in seq_len(ncol(prob))) {
    inc[, j, ] <- prob[, j] * .delay_cdf(
      family, eval_times, pwindow, loc[, j], shape_or_sd[, j], delay_max
    )
  }
  inc
}

# Daily hazards from cumulative incidence arrays evaluated on `eval_times`: the
# probability of each outcome during [t, t + 1) divided by the probability of
# still being at risk at t, which excludes every resolved case for the
# cause-specific hazard and only that outcome's cases for the subdistribution
# hazard.
.daily_hazards <- function(incidence, eval_times, times, type) {
  now <- match(times, eval_times)
  nxt <- match(times + 1, eval_times)
  resolved <- Reduce(`+`, incidence)
  lapply(incidence, function(inc) {
    at_risk <- if (type == "hazard") {
      1 - resolved[, , now, drop = FALSE]
    } else {
      1 - inc[, , now, drop = FALSE]
    }
    h <- (inc[, , nxt, drop = FALSE] - inc[, , now, drop = FALSE]) / at_risk
    h[at_risk <= 0] <- NA_real_
    h
  })
}
