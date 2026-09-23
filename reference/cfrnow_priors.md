# Priors for the mixture-cure CFR model

Priors are placed on the interpretable onset-to-death delay **mean** and
**standard deviation** (in days) and on the CFR, so they do not depend
on the chosen delay family. Defaults are the BDBV/Ebolavirus values used
by the Julia reference model: an onset-to-death delay centred on the
Isiro 2012 line-list reanalysis (mean ~= 12.75 d, sd ~= 7 d) and an
EVD/BDBV `Beta(6.6, 13.4)` CFR prior (mean ~= 0.33). The delay mean and
sd carry half-normal priors (truncated at zero).

## Usage

``` r
cfrnow_priors(
  delay_mean_mean = 12.75,
  delay_mean_sd = 3,
  delay_sd_mean = 7,
  delay_sd_sd = 2,
  cfr_a = 6.6,
  cfr_b = 13.4
)
```

## Arguments

- delay_mean_mean, delay_mean_sd:

  Half-normal prior on the delay mean (days).

- delay_sd_mean, delay_sd_sd:

  Half-normal prior on the delay sd (days).

- cfr_a, cfr_b:

  Beta prior on the case fatality ratio.

## Value

A named list of prior hyperparameters for [`fit_cfr()`](fit_cfr.md).

## Examples

``` r
cfrnow_priors(delay_mean_mean = 10)
#> $delay_mean_mean
#> [1] 10
#> 
#> $delay_mean_sd
#> [1] 3
#> 
#> $delay_sd_mean
#> [1] 7
#> 
#> $delay_sd_sd
#> [1] 2
#> 
#> $cfr_a
#> [1] 6.6
#> 
#> $cfr_b
#> [1] 13.4
#> 
```
