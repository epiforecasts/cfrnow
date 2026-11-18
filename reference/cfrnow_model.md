# Compile the mixture-cure CFR Stan model

Compiles `inst/stan/cfrnow.stan`, resolving the vendored primarycensored
functions via the package's Stan include path. The compiled model is
cached by cmdstanr, so repeated calls are cheap.

## Usage

``` r
cfrnow_model(...)
```

## Arguments

- ...:

  Passed to
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).

## Value

A `CmdStanModel` object.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- cfrnow_model()
} # }
```
