#' Real-time CFR as an `epidist` mixture-cure model
#'
#' `cfrnow` registers a mixture-cure survival model as an [epidist::epidist()]
#' model type. Each case is fatal with probability `cfr` and, when fatal, dies
#' at an onset-to-death delay; cases still alive at the observation cut-off are
#' right-censored, which corrects the downward bias of the naive deaths / cases
#' ratio in real time. Because the model is an `epidist` subclass, the CFR and
#' the delay both take `brms` formulas, e.g.
#' `epidist(data, bf(mu ~ 1, cfr ~ age), family = lognormal())`.
#'
#' With a `recovery_delay`, the fit becomes a two-outcome mixture-cure model
#' that also times recoveries: a non-fatal case recovers at a second delay, so a
#' recovered case contributes `(1 - cfr) f_R(r)` and an unresolved case
#' `cfr (1 - F_D(t)) + (1 - cfr)(1 - F_R(t))`. The recovery delay may use a
#' different family from the death delay.
#'
#' The delay distribution's location is `mu` (as `epidist` expects); the cure
#' probability `cfr` is an additional dpar with a logit link. Supported delay
#' families are `lognormal()`, `Gamma()` and `Weibull()`.
#'
#' @name cfrnow-cure-model
#' @keywords internal
NULL

# Outcome codes carried in the `outcome` vreal (match cfrnow's likelihood):
#   1 = observed death (timed)       -> log(cfr) + death delay density
#   2 = recorded recovery (timed)    -> log1m(cfr) + recovery delay density
#   3 = resolved non-death (untimed) -> log1m(cfr)  (cure factor)
#   0 = censored survivor            -> death-only: log1m(cfr F_D(t));
#                                       two-outcome: log_sum_exp mixture
.CURE_DEATH <- 1
.CURE_RECOVERY <- 2
.CURE_RESOLVED <- 3
.CURE_CENSORED <- 0

#' Build an `epidist_cure_model` from a line list or `prepare_cfr_data()` output
#'
#' Produces the one-row-per-case data frame the model needs: `y` (a delay for
#' deaths and recoveries, the follow-up time for survivors), an `outcome` code,
#' and the primary (`pwindow`) and secondary (`swindow`) censoring widths. When
#' the input carries timed recoveries the result is flagged for the two-outcome
#' fit via a `use_recovery` attribute.
#'
#' @param data Either the list returned by [prepare_cfr_data()] or a data frame
#'   already carrying `y`, `outcome`, `pwindow` and `swindow`.
#' @return A data frame of class `epidist_cure_model`.
#' @family cure_model
#' @export
as_epidist_cure_model <- function(data) {
  if (inherits(data, "cfrnow_data")) {
    data <- data$cases # one row per kept case
  }
  stopifnot(c("y", "outcome", "pwindow", "swindow") %in% names(data))
  data <- as.data.frame(data)
  class(data) <- c("epidist_cure_model", "data.frame")
  attr(data, "use_recovery") <- any(data$outcome == .CURE_RECOVERY)
  data
}

#' @method assert_epidist epidist_cure_model
#' @export
assert_epidist.epidist_cure_model <- function(data, ...) {
  cols <- c("y", "outcome", "pwindow", "swindow")
  codes <- c(.CURE_DEATH, .CURE_RECOVERY, .CURE_RESOLVED, .CURE_CENSORED)
  checkmate::assert_data_frame(data)
  checkmate::assert_names(names(data), must.include = cols)
  checkmate::assert_integerish(data$y)
  checkmate::assert_subset(data$outcome, codes)
  invisible(TRUE)
}

# Per-dpar links for a delay family: mu from the family, other (positive) params
# default to log when the base family carried no explicit link.
.delay_links <- function(family) {
  other <- setdiff(family$dpars, "mu")
  ol <- vapply(other, function(dp) {
    l <- family[[paste0("link_", dp)]]
    if (is.null(l) || !nzchar(l)) "log" else l
  }, character(1))
  c(family$link, ol)
}

.assert_delay_family <- function(family, what = "delay") {
  if (!family$family %in% c("lognormal", "gamma", "weibull")) {
    stop("cfrnow supports lognormal(), Gamma() and Weibull() ", what, "s only.",
      call. = FALSE
    )
  }
}

#' @method epidist_family_model epidist_cure_model
#' @export
epidist_family_model.epidist_cure_model <- function(data, family, ...) {
  .assert_delay_family(family)
  dpars <- c(family$dpars, "cfr")
  links <- c(.delay_links(family), "logit")
  if (isTRUE(attr(data, "use_recovery"))) {
    rfam <- attr(data, "recovery_family")
    if (is.null(rfam)) rfam <- family
    .assert_delay_family(rfam, "recovery delay")
    dpars <- c(dpars, paste0("r", rfam$dpars)) # recovery: own family
    links <- c(links, .delay_links(rfam))
  }
  brms::custom_family(
    paste0("cfrnow_", family$family),
    dpars = dpars, links = links, type = "int",
    vars = c("vreal1[n]", "vreal2[n]", "vreal3[n]", "primary_params"),
    loop = TRUE,
    log_lik = function(i, prep, ...) rep(NA_real_, prep$ndraws),
    posterior_predict = function(i, prep, ...) rep(NA_real_, prep$ndraws)
  )
}

#' @method epidist_formula_model epidist_cure_model
#' @export
epidist_formula_model.epidist_cure_model <- function(data, formula, ...) {
  stats::update(formula, y | vreal(outcome, pwindow, swindow) ~ .)
}

#' @method epidist_transform_data_model epidist_cure_model
#' @export
epidist_transform_data_model.epidist_cure_model <- function(data, family,
                                                            formula, ...) {
  data
}

#' @method epidist_model_prior epidist_cure_model
#' @export
epidist_model_prior.epidist_cure_model <- function(data, formula, ...) NULL

# The native Stan parameterisation brms uses for a family's mu-form, read out of
# a dummy model (as epidist does), so primarycensored gets the right arguments.
# The mu-form arguments may themselves contain parentheses (e.g. weibull's
# `mu ./ tgamma(1 + 1 ./ shape)`), so the args run to the statement's closing
# `);` rather than the first `)`.
.family_stan_param <- function(family_name) {
  code <- brms::make_stancode(y ~ 1,
    data = data.frame(y = c(1, 2)),
    family = family_name
  )
  pat <- sprintf("%s_l[a-z]+\\(Y \\| (.+?)\\);", family_name)
  stmts <- regmatches(code, gregexpr(pat, code, perl = TRUE))[[1]]
  stmt <- grep("mu", stmts, fixed = TRUE, value = TRUE)[1]
  sub(pat, "\\1", stmt, perl = TRUE)
}

# Fill a `<<hole>>` template from inst/stan with a named list of replacements.
.fill_stan_template <- function(file, holes) {
  path <- system.file("stan", file, package = "cfrnow")
  code <- paste(readLines(path), collapse = "\n")
  for (nm in names(holes)) {
    code <- gsub(paste0("<<", nm, ">>"), holes[[nm]], code, fixed = TRUE)
  }
  code
}

# Template holes for the recovery half of a two-outcome fit: its own family's
# parameter declarations, distribution id, and native reparameterisation
# (r-prefixed to match the recovery dpars rmu, rsigma / rshape).
.recovery_holes <- function(family, recovery_family, recovery_dpars) {
  rfam <- recovery_family %||% family
  reparam <- .family_stan_param(rfam$family)
  for (dp in rfam$dpars) {
    reparam <- gsub(paste0("\\b", dp, "\\b"), paste0("r", dp), reparam)
  }
  rid <- primarycensored::pcd_stan_dist_id(rfam$family, type = "delay")
  list(
    recovery_pars = toString(paste0("real ", recovery_dpars)),
    recovery_id = rid,
    recovery_reparam = reparam
  )
}

#' @method epidist_stancode epidist_cure_model
#' @export
epidist_stancode.epidist_cure_model <- function(data, family, formula, ...) {
  family_name <- sub("^cfrnow_", "", family$name)
  cfr_pos <- match("cfr", family$dpars)
  delay_dpars <- family$dpars[seq_len(cfr_pos - 1)]
  use_recovery <- cfr_pos < length(family$dpars)

  holes <- list(
    family = family_name,
    death_pars = toString(paste0("real ", delay_dpars)),
    death_id = primarycensored::pcd_stan_dist_id(family_name, type = "delay"),
    death_reparam = family$param,
    primary_id = primarycensored::pcd_stan_dist_id("uniform", type = "primary")
  )
  template <- "cure_lpmf_death.stan"
  if (use_recovery) {
    recovery_dpars <- family$dpars[(cfr_pos + 1):length(family$dpars)]
    holes <- c(holes, .recovery_holes(
      family, attr(data, "recovery_family"),
      recovery_dpars
    ))
    template <- "cure_lpmf_two_outcome.stan"
  }
  lpmf <- .fill_stan_template(template, holes)

  brms::stanvar(
    scode = primarycensored::pcd_load_stan_functions(),
    block = "functions"
  ) +
    brms::stanvar(scode = lpmf, block = "functions") +
    brms::stanvar(scode = "array[0] real primary_params;", block = "parameters")
}
