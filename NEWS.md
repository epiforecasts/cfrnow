# cfrnow (development version)

* `fit_cfr()` gains `loss_prior`, a `Beta()` prior on the probability that a
  case is lost to follow-up and its outcome never recorded. A fit that allows
  for loss estimates it alongside the outcome probability and the delays, and
  `summary()` reports it as a `loss` row. Without it, cases that stay
  unresolved far longer than the delays allow either stretch the delay's tail
  or stop the sampler from starting.
* `prepare_cfr_data()` gains `last_contact_date`, the column holding the date a
  case with no recorded outcome was last known unresolved. Such a case is
  censored there instead of at the cut-off, in retrospective fits too.
* The survival term for an unresolved case is computed on the log scale
  (`primarycensored_lcdf()` and `log1m_exp()`). A case followed up for much
  longer than the delay rounds the CDF to 1, where the previous `log1m()` form
  rejected every draw with `log1m: x is 1, but must be less than or equal to
  1`, leaving the chains stuck at their starting values.

* A delay's `max` is now honoured: `fit_cfr()` truncates the fitted delay at
  the bound (`LogNormal(..., max = 30)`), `simulate_linelist()` draws from the
  truncated delay, and `pp_check_cfr()` replicates from it. Previously the bound
  was silently ignored. The bound applies to the recorded delay, so it runs from
  0 to `max - 1` days, the same support distspec gives the delay object.
  `fit_cfr()` also stops, with a message naming the cases, when a recorded delay
  or an unresolved case falls outside the bounds and would otherwise fail inside
  Stan.

* `pp_check_cfr()` now draws replicate delays from a Weibull fit's own
  distribution; it previously drew them from a gamma, so the checks for a
  Weibull delay compared the fit against the wrong replicates.

* The model parameter (and its prior) is renamed from `cfr` to `prob`, since it
  is a case fatality ratio only when the line list runs from onset to death; the
  same model can fit a hospital fatality ratio or other outcome probability for a
  differently defined line list. Use `prob ~ ...` in `formula` and
  `prob_prior` in `fit_cfr()`; `summary()` now reports a `prob` (or
  `prob[<group>]`) row. `cfr ~ ...` and `cfr_prior` are still accepted and
  translated to `prob`, with a soft-deprecation warning.

# cfrnow 0.2.1

* `fit_cfr()` warns and `pp_check_cfr()` no longer errors when a `formula`
  covariate has missing values: brms drops those cases before fitting, and
  the stored onset dates are now kept in step with the rows it actually used.
* distspec is now on CRAN, so it is dropped from `Remotes` and installed from
  CRAN like the other dependencies.

# cfrnow 0.2.0

* `fit_cfr()` and `simulate_linelist()` support a `Weibull()` onset-to-death (and
  recovery) delay, alongside `LogNormal()` and `Gamma()`.
* Delay parameterisation now uses distspec's exported `natural_params()` in place
  of an internal helper, tracking the distspec API.
* Added a "Stratified and partially-pooled CFR" vignette covering no-, complete-
  and partial-pooling CFR fits and per-group `summary()` output.
* `pp_check_cfr()` runs a posterior-predictive check on a fit: it draws replicate
  line-list outcomes from the posterior, replays the real-time truncation, and
  compares the observed death counts (plus recoveries in a two-outcome fit) and
  the observed onset-to-death delays against the replicates (#14).
* `summary()` gains an `ascertainment_ratio` argument that corrects the CFR for
  outcome-dependent case ascertainment (fatal and non-fatal cases entering the
  line list at different rates). The ratio is supplied, defaulting to 1.
* `fit_cfr()` accepts intercept-free CFR formulas (e.g. `cfr ~ 0 + group`, one
  estimated logit-CFR per group): the `cfr_prior` is placed on those
  coefficients rather than a non-existent intercept, so the fit no longer fails
  brms prior validation.
* `summary()` reports a CFR per group for a `cfr ~ group` fit (one `cfr[<group>]`
  row per group), rather than erroring or silently reporting only the reference
  level.

# cfrnow 0.1.0

First release.

* `fit_cfr()` estimates a real-time case fatality ratio from line-list data with a
  Bayesian mixture-cure survival model. It is registered as an `epidist` model
  type, so the CFR and the onset-to-death delay both take `brms` formulas.
* `prepare_cfr_data()` turns a line list into model inputs. It sorts each case, at
  a chosen observation cut-off, into an observed death, a resolved non-death, or a
  right-censored survivor.
* The onset-to-death delay (LogNormal or Gamma) can be co-estimated or held fixed.
  Hold it fixed and you get the Ghani/Nishiura estimator.
* Pass a `recovery_date` column and a two-outcome fit also times recoveries.
* Put a `brms` formula on the CFR or the delay for covariates or a time-varying
  CFR.
* `simulate_linelist()` builds line lists for testing and examples.
* `summary()` and `print()` report the corrected CFR, the delay moments,
  convergence diagnostics, and a flag for when the CFR is only weakly identified.
