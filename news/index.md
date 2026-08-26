# Changelog

## cfrnow 0.2.1

- distspec is now on CRAN, so it is dropped from `Remotes` and installed
  from CRAN like the other dependencies.

## cfrnow 0.2.0

- [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md) and
  [`simulate_linelist()`](https://epiforecasts.io/cfrnow/reference/simulate_linelist.md)
  support a
  [`Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html)
  onset-to-death (and recovery) delay, alongside
  [`LogNormal()`](https://epiforecasts.io/distspec/reference/LogNormal.html)
  and
  [`Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html).
- Delay parameterisation now uses distspec’s exported
  [`natural_params()`](https://epiforecasts.io/distspec/reference/natural_params.html)
  in place of an internal helper, tracking the distspec API.
- Added a “Stratified and partially-pooled CFR” vignette covering no-,
  complete- and partial-pooling CFR fits and per-group
  [`summary()`](https://rdrr.io/r/base/summary.html) output.
- [`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md)
  runs a posterior-predictive check on a fit: it draws replicate
  line-list outcomes from the posterior, replays the real-time
  truncation, and compares the observed death counts (plus recoveries in
  a two-outcome fit) and the observed onset-to-death delays against the
  replicates ([\#14](https://github.com/epiforecasts/cfrnow/issues/14)).
- [`summary()`](https://rdrr.io/r/base/summary.html) gains an
  `ascertainment_ratio` argument that corrects the CFR for
  outcome-dependent case ascertainment (fatal and non-fatal cases
  entering the line list at different rates). The ratio is supplied,
  defaulting to 1.
- [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md)
  accepts intercept-free CFR formulas (e.g. `cfr ~ 0 + group`, one
  estimated logit-CFR per group): the `cfr_prior` is placed on those
  coefficients rather than a non-existent intercept, so the fit no
  longer fails brms prior validation.
- [`summary()`](https://rdrr.io/r/base/summary.html) reports a CFR per
  group for a `cfr ~ group` fit (one `cfr[<group>]` row per group),
  rather than erroring or silently reporting only the reference level.

## cfrnow 0.1.0

First release.

- [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md)
  estimates a real-time case fatality ratio from line-list data with a
  Bayesian mixture-cure survival model. It is registered as an `epidist`
  model type, so the CFR and the onset-to-death delay both take `brms`
  formulas.
- [`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
  turns a line list into model inputs. It sorts each case, at a chosen
  observation cut-off, into an observed death, a resolved non-death, or
  a right-censored survivor.
- The onset-to-death delay (LogNormal or Gamma) can be co-estimated or
  held fixed. Hold it fixed and you get the Ghani/Nishiura estimator.
- Pass a `recovery_date` column and a two-outcome fit also times
  recoveries.
- Put a `brms` formula on the CFR or the delay for covariates or a
  time-varying CFR.
- [`simulate_linelist()`](https://epiforecasts.io/cfrnow/reference/simulate_linelist.md)
  builds line lists for testing and examples.
- [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`print()`](https://rdrr.io/r/base/print.html) report the corrected
  CFR, the delay moments, convergence diagnostics, and a flag for when
  the CFR is only weakly identified.
