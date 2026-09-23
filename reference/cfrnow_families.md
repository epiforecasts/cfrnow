# Supported onset-to-death delay families

The families `cfrnow` can fit, each an analytical primarycensored delay
distribution paired with a uniform primary event. These match the
dist.spec constructors accepted by
[`fit_cfr()`](https://sbfnk.github.io/cfrnow/reference/fit_cfr.md)'s
`delay` argument
([`dist.spec::LogNormal()`](https://epiforecasts.io/dist.spec/reference/Distributions.html),
[`dist.spec::Gamma()`](https://epiforecasts.io/dist.spec/reference/Distributions.html)).

## Usage

``` r
cfrnow_families()
```

## Value

A character vector of supported family names.

## Examples

``` r
cfrnow_families()
#> [1] "lognormal" "gamma"    
```
