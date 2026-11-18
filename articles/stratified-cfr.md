# Stratified and partially-pooled CFR

## Why stratify, and why pool

The CFR often varies by group – age, treatment centre, calendar period –
and a per-group estimate is frequently required. Two limiting cases
bracket the options. Fitting each group independently (*no pooling*)
gives high-variance estimates where a group’s sample is small. Fitting
one CFR for all groups (*complete pooling*) ignores between-group
variation.

*Partial pooling* lies between the two: each group’s CFR is drawn from a
common distribution, so groups with more data are estimated close to
independently while groups with less data are shrunk toward the shared
mean. `cfrnow` is an [epidist](https://epidist.epinowcast.org/) model,
so `prob` takes a `brms` formula and pooling is a one-line change:
`prob ~ (1 | group)`.

## A six-site line list

We simulate six treatment centres: two large (n = 450, 250) and four
small (n = 15–20), with true CFRs spread around a common mean. The large
sites pin the shared mean; the small ones carry little information on
their own, so shrinkage acts mainly on them. Site labels are synthetic.

``` r

library(cfrnow)
#> Loading required package: distspec
#> 
#> Attaching package: 'distspec'
#> The following objects are masked from 'package:stats':
#> 
#>     Gamma, sd
set.seed(20260714)

sites <- data.frame(
  site = c("A", "B", "C", "D", "E", "F"),
  n    = c(450, 250, 20, 15, 15, 18),
  cfr  = c(0.45, 0.50, 0.30, 0.65, 0.35, 0.60)
)

death_delay    <- LogNormal(mean = 11, sd = 6.5)
recovery_delay <- LogNormal(mean = 16, sd = 8)

ll <- do.call(rbind, Map(function(s, n, cfr) {
  cbind(
    simulate_linelist(
      n = n, cfr = cfr, delay = death_delay,
      recovery = recovery_delay, onset_days = 45
    ),
    site = s
  )
}, sites$site, sites$n, sites$cfr))

d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-20"), covariates = "site")
c(cases = d$n_cases, deaths = d$n_deaths)
#>  cases deaths 
#>    768    331
```

## Three fits: complete, no, and partial pooling

The only thing that changes across the three fits is the CFR formula;
the delay and recovery delay are co-estimated in each. (Two short chains
keep the build time short; use more for inference.)

``` r

otd <- LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15))
otr <- LogNormal(meanlog = Normal(2.7, 0.2), sdlog = Normal(0.5, 0.15))
args <- list(
  delay = otd, recovery_delay = otr, prob_prior = Beta(1, 1),
  backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 1
)

# complete pooling: a single CFR
fit_pool    <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, prob ~ 1)), args))

# no pooling: an independent CFR per site (one estimated logit-CFR per site)
fit_sep     <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, prob ~ 0 + site)), args))

# partial pooling: per-site CFR shrunk to a common mean
fit_partial <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, prob ~ (1 | site))), args))
```

## Reading the estimates

For the complete-pool fit,
[`summary()`](https://rdrr.io/r/base/summary.html) reports a single
`prob` row:

``` r

summary(fit_pool)
#>        quantity      mean       q2.5        q50      q97.5      rhat  ess_bulk
#> 1          prob  0.506545  0.4636737  0.5071414  0.5485087 0.9983105  830.2913
#> 2    delay_mean 11.787311 10.7885060 11.7894961 12.8923528 1.0029524 1056.9733
#> 3      delay_sd  8.221317  6.8659083  8.2153459  9.7142108 1.0035616 1041.1778
#> 4 recovery_mean 15.417467 14.4281222 15.4091710 16.5072024 1.0008222  749.2406
#> 5   recovery_sd  7.780150  6.7430025  7.7225057  9.1285643 1.0010545  749.5455
```

For a grouped fit it reports one `prob[<site>]` row per site:

``` r

summary(fit_sep)
#>         quantity       mean        q2.5        q50      q97.5      rhat
#> 1        prob[A]  0.4791414  0.42061391  0.4799069  0.5335332 1.0024523
#> 2        prob[B]  0.5638484  0.49978464  0.5644529  0.6305058 1.0006084
#> 3        prob[C]  0.1915139  0.05379895  0.1825094  0.3822915 1.0003174
#> 4        prob[D]  0.5191643  0.27106838  0.5217161  0.7575433 1.0009786
#> 5        prob[E]  0.5910984  0.29094955  0.5978461  0.8619616 1.0057751
#> 6        prob[F]  0.7712562  0.51710526  0.7818849  0.9511256 1.0040599
#> 7     delay_mean 11.8565652 10.89093661 11.8394599 12.9353615 0.9995752
#> 8       delay_sd  8.3060791  7.00723709  8.2798773  9.8823679 0.9990411
#> 9  recovery_mean 15.3554429 14.43829181 15.3235797 16.4154553 1.0011378
#> 10   recovery_sd  7.7211795  6.63175154  7.6981082  8.8849930 1.0053201
#>    ess_bulk
#> 1  1348.233
#> 2  1993.419
#> 3  2219.770
#> 4  2408.730
#> 5  1921.446
#> 6  1824.690
#> 7  1340.579
#> 8  1333.187
#> 9  1429.659
#> 10 1257.175
```

Collect those per-site CFR rows from the two grouped fits for the figure
below; the complete-pool fit supplies the shared mean.

``` r

cfr_rows <- function(fit, approach) {
  s <- summary(fit)
  r <- s[grepl("^prob\\[", s$quantity), ] # the prob[<site>] rows
  data.frame(
    approach = approach,
    site = sub("^prob\\[(.*)\\]$", "\\1", r$quantity),
    median = r$q50, lower = r$q2.5, upper = r$q97.5
  )
}

est <- rbind(cfr_rows(fit_sep, "no pool"), cfr_rows(fit_partial, "partial pool"))
est <- merge(est, sites[, c("site", "n", "cfr")], by = "site")
names(est)[names(est) == "cfr"] <- "truth"

s_pool <- summary(fit_pool)
pool_mean <- s_pool$q50[s_pool$quantity == "prob"]
```

## Shrinkage

Each site shows its no-pool and partial-pool estimate against the truth
(×) and the complete-pool mean (dotted). For the large sites (A, B) the
two estimates nearly coincide. For the small sites (C–F) the no-pool
interval is wide and the estimate noisy; partial pooling pulls it toward
the mean and narrows the interval sharply – lower variance at the cost
of some bias.

``` r

library(ggplot2)
est$approach <- factor(est$approach, levels = c("no pool", "partial pool"))
est$site <- factor(est$site, levels = rev(sites$site))

ggplot(est, aes(median, site, colour = approach)) +
  geom_vline(xintercept = pool_mean, linetype = 3) +
  geom_point(aes(x = truth), shape = 4, size = 2.4, colour = "grey30") +
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.55)) +
  scale_y_discrete(labels = function(s) paste0(s, "  (n=", sites$n[match(s, sites$site)], ")")) +
  labs(
    x = "CFR: median & 95% CrI   (× = truth, dotted = complete-pool mean)",
    y = NULL, colour = NULL
  ) +
  theme_minimal()
```

![](stratified-cfr_files/figure-html/shrinkage-1.png)

## Varying the delay by group

A pooled onset-to-death delay that appears irregular is often a mixture
of several populations. If sites resolve on different timescales, pool
the delay by site as well, so each retains its own onset-to-death
distribution rather than contributing to a single multimodal one:

``` r

fit_both <- do.call(fit_cfr, c(
  list(d, formula = brms::bf(mu ~ (1 | site), prob ~ (1 | site))), args
))
```

## Check the fit

Check that the parametric delay reproduces the observed delays; the
real-time correction depends on extrapolating that delay’s tail.

``` r

pp_check_cfr(fit_partial, type = "delay")
```

![](stratified-cfr_files/figure-html/ppcheck-1.png)
