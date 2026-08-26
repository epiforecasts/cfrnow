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
so the CFR takes a `brms` formula and pooling is a one-line change:
`cfr ~ (1 | group)`.

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
#>    768    313
```

## Three fits: complete, no, and partial pooling

The only thing that changes across the three fits is the CFR formula;
the delay and recovery delay are co-estimated in each. (Two short chains
keep the build time short; use more for inference.)

``` r

otd <- LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15))
otr <- LogNormal(meanlog = Normal(2.7, 0.2), sdlog = Normal(0.5, 0.15))
args <- list(
  delay = otd, recovery_delay = otr, cfr_prior = Beta(1, 1),
  backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 1
)

# complete pooling: a single CFR
fit_pool    <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, cfr ~ 1)), args))

# no pooling: an independent CFR per site (one estimated logit-CFR per site)
fit_sep     <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, cfr ~ 0 + site)), args))

# partial pooling: per-site CFR shrunk to a common mean
fit_partial <- do.call(fit_cfr, c(list(d, formula = brms::bf(mu ~ 1, cfr ~ (1 | site))), args))
```

## Reading the estimates

For the complete-pool fit,
[`summary()`](https://rdrr.io/r/base/summary.html) reports a single
`cfr` row:

``` r

summary(fit_pool)
#>        quantity       mean       q2.5        q50      q97.5      rhat  ess_bulk
#> 1           cfr  0.4678748  0.4286629  0.4673893  0.5066603 0.9990148 1164.7840
#> 2    delay_mean 10.8730852 10.1083283 10.8453649 11.8350403 1.0004768  962.9530
#> 3      delay_sd  6.7168988  5.6859892  6.6600172  8.0109392 0.9996777 1021.0216
#> 4 recovery_mean 16.0704570 15.1006704 16.0484178 17.1530976 1.0008807  867.6661
#> 5   recovery_sd  8.2813389  7.1671840  8.2444780  9.5517492 1.0019728  827.1322
```

For a grouped fit it reports one `cfr[<site>]` row per site:

``` r

summary(fit_sep)
#>         quantity       mean       q2.5        q50      q97.5      rhat ess_bulk
#> 1         cfr[A]  0.4474581  0.3960099  0.4480130  0.4992457 0.9985338 2470.609
#> 2         cfr[B]  0.4863861  0.4175876  0.4855172  0.5583176 1.0000843 1777.778
#> 3         cfr[C]  0.4265776  0.2165792  0.4266473  0.6469009 1.0134665 1707.025
#> 4         cfr[D]  0.6816893  0.4295505  0.6907186  0.8891901 0.9993289 1924.529
#> 5         cfr[E]  0.5141502  0.2614512  0.5141866  0.7649760 1.0019766 1441.862
#> 6         cfr[F]  0.5369987  0.2872178  0.5416325  0.7707196 1.0018319 2741.035
#> 7     delay_mean 10.8779669 10.0646285 10.8728879 11.8695340 0.9996002 1722.205
#> 8       delay_sd  6.7248734  5.7165388  6.6646874  8.0242228 0.9998720 1773.031
#> 9  recovery_mean 16.0588590 15.0737446 16.0279880 17.0938482 1.0000791 1789.818
#> 10   recovery_sd  8.2676158  7.2480353  8.2257169  9.4818801 1.0031554 1820.836
```

Collect those per-site CFR rows from the two grouped fits for the figure
below; the complete-pool fit supplies the shared mean.

``` r

cfr_rows <- function(fit, approach) {
  s <- summary(fit)
  r <- s[grepl("^cfr\\[", s$quantity), ] # the cfr[<site>] rows
  data.frame(
    approach = approach,
    site = sub("^cfr\\[(.*)\\]$", "\\1", r$quantity),
    median = r$q50, lower = r$q2.5, upper = r$q97.5
  )
}

est <- rbind(cfr_rows(fit_sep, "no pool"), cfr_rows(fit_partial, "partial pool"))
est <- merge(est, sites[, c("site", "n", "cfr")], by = "site")
names(est)[names(est) == "cfr"] <- "truth"

s_pool <- summary(fit_pool)
pool_mean <- s_pool$q50[s_pool$quantity == "cfr"]
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
  list(d, formula = brms::bf(mu ~ (1 | site), cfr ~ (1 | site))), args
))
```

## Check the fit

Check that the parametric delay reproduces the observed delays; the
real-time correction depends on extrapolating that delay’s tail.

``` r

pp_check_cfr(fit_partial, type = "delay")
```

![](stratified-cfr_files/figure-html/ppcheck-1.png)
