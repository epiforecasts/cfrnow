# Real-time CFR as an `epidist` mixture-cure model

`cfrnow` registers a mixture-cure survival model as an
[`epidist::epidist()`](https://epidist.epinowcast.org/reference/epidist.html)
model type. Each case is fatal with probability `prob` and, when fatal,
dies at an onset-to-death delay; cases still alive at the observation
cut-off are right-censored, which corrects the downward bias of the
naive deaths / cases ratio in real time. Because the model is an
`epidist` subclass, `prob` and the delay both take `brms` formulas, e.g.
`epidist(data, bf(mu ~ 1, prob ~ age), family = lognormal())`. `prob` is
a CFR when the line list runs from onset to death; the same model fits a
hospital fatality ratio or other outcome probability for a differently
defined line list.

## Details

With a `recovery_delay`, the fit becomes a two-outcome mixture-cure
model that also times recoveries: a non-fatal case recovers at a second
delay, so a recovered case contributes `(1 - prob) f_R(r)` and an
unresolved case `prob (1 - F_D(t)) + (1 - prob)(1 - F_R(t))`. The
recovery delay may use a different family from the death delay.

The delay distribution's location is `mu` (as `epidist` expects); the
cure probability `prob` is an additional dpar with a logit link.
Supported delay families are `lognormal()`,
[`Gamma()`](https://epiforecasts.io/distspec/reference/Gamma.html) and
[`Weibull()`](https://epiforecasts.io/distspec/reference/Weibull.html).
