# curecfr: Real-Time Case Fatality Ratio via a Mixture-Cure Survival Model

Estimates a real-time case fatality ratio (CFR) from line-list data
using a Bayesian mixture-cure survival model. Each case is fatal with
probability equal to the CFR and, when fatal, dies at an
interval-censored onset-to-death delay; cases still alive at the
observation cut-off are right-censored, which corrects the downward bias
of the naive deaths / cases ratio in real time. The onset-to-death delay
is co-estimated using the primarycensored analytical
censored-distribution machinery.

## Author

**Maintainer**: Sebastian Funk <sebastian.funk@lshtm.ac.uk>

Authors:

- Sebastian Funk <sebastian.funk@lshtm.ac.uk>
