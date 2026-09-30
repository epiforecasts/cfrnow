# Default onset-to-death delay specification

A
[`dist.spec::LogNormal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
with
[`Normal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
priors on the native parameters, reproducing the BDBV/Ebolavirus
onset-to-death prior of the Julia reference model (the Isiro 2012
line-list reanalysis): `meanlog ~ Normal(2.41, 0.2)`,
`sdlog ~ Normal(0.51, 0.15)` (delay mean ~= 12.75 d, sd ~= 7 d).

## Usage

``` r
default_delay()
```

## Value

A `dist_spec` for the onset-to-death delay.

## Details

Supply your own
[`dist.spec::LogNormal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
or
[`dist.spec::Gamma()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
to [`fit_cfr()`](https://sbfnk.github.io/cfrnow/reference/fit_cfr.md) to
change the family or priors. Give each native parameter as a
[`Normal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
to co-estimate it, or as a fixed number /
[`dist.spec::Fixed()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)
to hold it fixed (fixing the whole delay gives the Ghani/Nishiura
estimator).

## Examples

``` r
default_delay()
#> - lognormal distribution:
#>   meanlog:
#>     - normal distribution:
#>       mean:
#>         2.4
#>       sd:
#>         0.2
#>   sdlog:
#>     - normal distribution:
#>       mean:
#>         0.51
#>       sd:
#>         0.15
```
