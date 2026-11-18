# Posterior cumulative incidence and hazards of death and recovery

Turns a fit's outcome probability and delay distributions into the
quantities a competing-risks analysis reports, for each covariate
pattern in `newdata`. The cumulative incidence of death and of recovery
by time `t` since onset are \$\$I_D(t) = p \\ F_D(t), \qquad I_R(t) =
(1 - p) \\ F_R(t)\$\$ where `p` is `prob` and `F_D` and `F_R` are the
distribution functions of the death and recovery delays, censored for
the day of onset and truncated at each delay's `max` as in the
likelihood.

## Usage

``` r
posterior_cumulative_incidence(
  object,
  newdata = NULL,
  times = NULL,
  type = c("cumulative_incidence", "hazard", "subdistribution_hazard"),
  pwindow = 1
)
```

## Arguments

- object:

  A `cfrnow_fit` from
  [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md).

- newdata:

  A data frame with the covariates of the `prob` and delay formulas, one
  row per covariate pattern. Defaults to the distinct covariate patterns
  among the fitted cases, in the order they first appear in
  `object$data` (a single row for a fit without covariates).

- times:

  Days since onset at which to evaluate. Defaults to every day from 0 to
  the longest follow-up in the fitted data.

- type:

  `"cumulative_incidence"` (the default), `"hazard"` for the
  cause-specific hazard or `"subdistribution_hazard"`.

- pwindow:

  Width of the onset window, in days, that the delays are censored over.
  Defaults to a single day.

## Value

A data frame with one row per posterior draw, row of `newdata`, time and
outcome: `.draw`, `row` (the row of `newdata`), `time`, `outcome`
(`"death"` or `"recovery"`) and `value`. A hazard is `NA` at a time when
every case has already resolved.

## Details

The hazards are daily, matching the daily resolution of the data. The
cause-specific hazard of death on day `t` is the probability of dying
during `[t, t + 1)` among cases still unresolved at `t`, \$\$h_D(t) =
\frac{I_D(t + 1) - I_D(t)}{1 - I_D(t) - I_R(t)},\$\$ and the
subdistribution (Fine-Gray) hazard keeps cases that have recovered in
the risk set, \$\$\tilde h_D(t) = \frac{I_D(t + 1) - I_D(t)}{1 -
I_D(t)}.\$\$ Recovery works the same way. Hazard ratios between two
covariate patterns are ratios of these draws. They vary with `t`, since
the model does not assume proportional hazards; covariates act on `prob`
and on the delays separately.

A fit without a `recovery_delay` returns death only. Its survivors never
resolve in the model, so the cause-specific and subdistribution hazards
coincide. Loss to follow-up only changes which outcomes get recorded and
plays no part here.

## See also

Other fit:
[`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md),
[`posterior_prob_death()`](https://epiforecasts.io/cfrnow/reference/posterior_prob_death.md),
[`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md),
[`summary.cfrnow_fit()`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
ll$group <- sample(c("a", "b"), nrow(ll), replace = TRUE)
d <- prepare_cfr_data(ll,
  obs_time = max(ll$onset_date) - 5,
  covariates = "group"
)
fit <- fit_cfr(d,
  delay = LogNormal(Normal(2.4, 0.2), Normal(0.5, 0.15)),
  formula = brms::bf(mu ~ group, prob ~ group)
)
groups <- data.frame(group = c("a", "b"))

# cumulative incidence of death by day 7, 14 and 28 in each group
ci <- posterior_cumulative_incidence(fit, groups, times = c(7, 14, 28))
aggregate(value ~ row + time, ci, stats::median)

# cause-specific hazard ratio of death, group b vs a, over time
h <- posterior_cumulative_incidence(fit, groups,
  times = 0:28,
  type = "hazard"
)
a <- h[h$outcome == "death" & h$row == 1, ]
b <- h[h$outcome == "death" & h$row == 2, ]
a$hr <- b$value / a$value
aggregate(hr ~ time, a, stats::median)
} # }
```
