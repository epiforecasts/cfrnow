# Real-time CFR with curecfr

## The problem

The naive case fatality ratio, `deaths / cases`, is biased **downward**
in real time: recent cases have not yet had time to die. Restricting to
cases with a *resolved* outcome, `deaths / (deaths + recoveries)`, swaps
that for an **upward** bias, because deaths resolve faster than
recoveries, so at any mid-outbreak cut-off the resolved set is enriched
for deaths.

`curecfr` sidesteps both by never conditioning on resolution. Each case
is fatal with probability `cfr`; a fatal case dies at an
interval-censored onset-to-death delay `F`; every case still alive at
the cut-off is right-censored, contributing the mixture-cure survival
term `1 - cfr * F(t)`. Recovery times are not used, so the resolved-set
enrichment cannot arise.

``` r

library(curecfr)
```

## A line list

Your line list needs an `onset_date` and a `death_date` (`NA` for cases
that have not died); an optional `onset_lower`/`onset_upper` window
widens the onset censoring. Here we simulate one with a known CFR of
0.55 sampled part-way through a growing epidemic.

``` r

set.seed(1)
ll <- simulate_linelist(n = 500, cfr = 0.55, delay_mean = 12.75, delay_sd = 7,
                        delay_family = "gamma", onset_days = 45)
head(ll)
#>   onset_date death_date
#> 1 2026-01-04 2026-01-13
#> 2 2026-02-08 2026-03-15
#> 3 2026-01-01       <NA>
#> 4 2026-02-03 2026-02-05
#> 5 2026-01-23       <NA>
#> 6 2026-02-12       <NA>
```

## Preparing the data at a cut-off

[`prepare_cfr_data()`](../reference/prepare_cfr_data.md) classifies each
case at the observation cut-off into an observed death, a right-censored
survivor, or (retrospectively) a fully-resolved non-death.

``` r

cutoff <- max(ll$onset_date) - 3
d <- prepare_cfr_data(ll, obs_time = cutoff)
#> 35 case(s) with onset after the cut-off excluded
c(cases = d$n_cases, deaths = d$n_deaths, censored = d$n_cens)
#>    cases   deaths censored 
#>      465      178      287
naive <- d$n_deaths / d$n_cases
round(naive, 3)  # the downward-biased naive ratio
#> [1] 0.383
```

## Fitting

[`fit_cfr()`](../reference/fit_cfr.md) fits the mixture-cure model with
cmdstanr. The delay is parameterised by its mean and standard deviation,
so the choice of family (`"gamma"` by default, or `"lognormal"`) leaves
the priors and summaries unchanged. This step needs a working CmdStan
install and is not run when the vignette is built.

``` r

fit <- fit_cfr(d, delay_family = "gamma")
summarise_cfr(fit)
#>     quantity   mean   q2.5    q50  q97.5  rhat ess_bulk
#> 1        cfr  0.579  0.528  0.579  0.625 1.001     3200
#> 2 delay_mean 12.700 11.500 12.700 13.900 1.000     3400
#> 3   delay_sd  7.200  6.300  7.200  8.200 1.000     3100
```

The corrected CFR (~0.58) recovers the true 0.55 and sits well above the
naive 0.38, which is depressed by the many recent, not-yet-resolved
cases. The summary also carries `rhat`/`ess_bulk` and a
`cfr_low_information` attribute (see below). Pass `obs_time = NULL` to
[`prepare_cfr_data()`](../reference/prepare_cfr_data.md) for a
retrospective, fully-resolved fit, which returns the naive ratio with
the delay as a nuisance.

## Assumptions and caveats

The correction fixes a *timing* bias; it cannot repair the data.

- **Complete death ascertainment.** A death missing from the line list
  (e.g. a community death outside a treatment centre) is treated as a
  survivor and biases the CFR down — the dominant risk when
  ascertainment is ETC-centred.
- **Deaths known on the day they occur.** Notification lag reintroduces
  the bias the model removes, via the data rather than the delay;
  nowcast or caveat the death series if that lag is non-trivial.
- **Weak identifiability early on.** When few deaths have resolved, the
  CFR posterior stays close to its prior;
  [`summarise_cfr()`](../reference/summarise_cfr.md) flags this via
  `cfr_low_information`. Treat a flagged estimate as prior-driven.
- **Stationary delay and homogeneous CFR.** A pooled CFR lags reality as
  treatment access scales up; time-varying CFR is future work.
- **Delay family is chosen, not estimated.** Refit with
  `delay_family = "lognormal"` to check sensitivity to the tail.
