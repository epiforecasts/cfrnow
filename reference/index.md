# Package index

## Prepare and fit

- [`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
  : Build model inputs from a line list

- [`as_epidist_cure_model()`](https://epiforecasts.io/cfrnow/reference/as_epidist_cure_model.md)
  :

  Build an `epidist_cure_model` from a line list or
  [`prepare_cfr_data()`](https://epiforecasts.io/cfrnow/reference/prepare_cfr_data.md)
  output

- [`fit_cfr()`](https://epiforecasts.io/cfrnow/reference/fit_cfr.md) :
  Fit the real-time mixture-cure CFR model

## Summarise and check

- [`summary(`*`<cfrnow_fit>`*`)`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)
  [`print(`*`<cfrnow_fit>`*`)`](https://epiforecasts.io/cfrnow/reference/summary.cfrnow_fit.md)
  : Summarise a mixture-cure CFR fit
- [`pp_check_cfr()`](https://epiforecasts.io/cfrnow/reference/pp_check_cfr.md)
  : Posterior-predictive check for a cfrnow fit

## Simulation

- [`simulate_linelist()`](https://epiforecasts.io/cfrnow/reference/simulate_linelist.md)
  : Simulate a line list for testing and examples
