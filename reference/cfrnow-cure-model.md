# Real-time CFR as an `epidist` mixture-cure model

`cfrnow` registers a mixture-cure survival model as an
[`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html)
model type. Each case is fatal with probability `cfr` and, when fatal,
dies at an onset-to-death delay; cases still alive at the observation
cut-off are right-censored, which corrects the downward bias of the
naive deaths / cases ratio in real time. Because the model is an
`epidist` subclass, the CFR and the delay both take `brms` formulas,
e.g. `epidist(data, bf(mu ~ 1, cfr ~ age), family = lognormal())`.

## Details

With a `recovery_delay`, the fit becomes a two-outcome mixture-cure
model that also times recoveries: a non-fatal case recovers at a second
delay, so a recovered case contributes `(1 - cfr) f_R(r)` and an
unresolved case `cfr (1 - F_D(t)) + (1 - cfr)(1 - F_R(t))`. The recovery
delay may use a different family from the death delay.

The delay distribution's location is `mu` (as `epidist` expects); the
cure probability `cfr` is an additional dpar with a logit link.
Supported delay families are `lognormal()`,
[`Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html) and
[`Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html).
