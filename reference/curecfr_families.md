# Supported onset-to-death delay families

The families `curecfr` can fit, each an analytical primarycensored delay
distribution paired with a uniform primary event and parameterised here
by its mean and standard deviation. Weibull is analytically available in
primarycensored but not yet wired into the mean/sd conversion, so it is
excluded here.

## Usage

``` r
curecfr_families()
```

## Value

A character vector of supported family names.

## Examples

``` r
curecfr_families()
#> [1] "lognormal" "gamma"    
```
