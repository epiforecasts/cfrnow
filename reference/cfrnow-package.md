# cfrnow: Real-Time Case Fatality Ratio via a Mixture-Cure Survival Model

Estimates a real-time case fatality ratio (CFR) from line-list data
using a Bayesian mixture-cure survival model, registered as an 'epidist'
model type so the CFR and the onset-to-death delay both take 'brms'
formulas. Each case is fatal with probability equal to the CFR and, when
fatal, dies at an interval-censored onset-to-death delay; cases still
alive at the observation cut-off are right-censored, which corrects the
downward bias of the naive deaths / cases ratio in real time. The delay
is co-estimated using the 'primarycensored' censored-distribution
machinery.

## See also

Useful links:

- <https://epiforecasts.io/cfrnow/>

- <https://github.com/epiforecasts/cfrnow>

- Report bugs at <https://github.com/epiforecasts/cfrnow/issues>

## Author

**Maintainer**: Sebastian Funk <sebastian.funk@lshtm.ac.uk>
