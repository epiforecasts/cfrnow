# Default onset-to-recovery delay specification

A weakly-informative
[`dist.spec::LogNormal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
for the onset-to-recovery (discharge) delay, centred around a mean of
~20 days and sd of ~11 days, longer than onset-to-death. Pass it (or
your own) as `recovery_delay` to
[`fit_cfr()`](https://sbfnk.github.io/cfrnow/reference/fit_cfr.md) to
switch on the competing-risks model that uses recovery *timing*. There
is much less external information on this delay than on onset-to-death,
so set it from your own data where possible.

## Usage

``` r
default_recovery_delay()
```

## Value

A `dist_spec` for the onset-to-recovery delay.

## Examples

``` r
default_recovery_delay()
#> - lognormal distribution:
#>   meanlog:
#>     - normal distribution:
#>       mean:
#>         2.9
#>       sd:
#>         0.3
#>   sdlog:
#>     - normal distribution:
#>       mean:
#>         0.5
#>       sd:
#>         0.2
```
