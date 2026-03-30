# Violin plots of fixed effect coefficients in a compnet model

Violin plots of fixed effect coefficients in a compnet model

## Usage

``` r
fixedeff_violins(mod)
```

## Arguments

- mod:

  An object of class "compnet" created by the buildcompnet() function.

## Value

A ggplot2 graphic showing violin plots of standarized effect sizes for
fixed effects in a compnet model.

## Examples

``` r
data(ex_presabs)
data(ex_traits)

# Quick demo run. Will prompt warnings.
# Run with default warmup and iter for good posterior sampling.
ex_compnet <- buildcompnet(presabs=ex_presabs,
spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#> [1] "You are currently running a compnet model with a Fisher's noncentral hypergeometric likelihood. This is the default option because there is strong theory supporting it. However, choosing a binomial likelihood (i.e., setting family='binomial') instead may result in a substantially faster run. This alternative option performs equally well in simulation-based testing. See https://kyle-rosenblad.github.io/compnet/ for more details"
#> 
#> SAMPLING FOR MODEL 'srm_fnchypg' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000382 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.82 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration:  1 / 20 [  5%]  (Warmup)
#> Chain 1: Iteration:  2 / 20 [ 10%]  (Warmup)
#> Chain 1: Iteration:  4 / 20 [ 20%]  (Warmup)
#> Chain 1: Iteration:  6 / 20 [ 30%]  (Warmup)
#> Chain 1: Iteration:  8 / 20 [ 40%]  (Warmup)
#> Chain 1: Iteration: 10 / 20 [ 50%]  (Warmup)
#> Chain 1: Iteration: 11 / 20 [ 55%]  (Sampling)
#> Chain 1: Iteration: 12 / 20 [ 60%]  (Sampling)
#> Chain 1: Iteration: 14 / 20 [ 70%]  (Sampling)
#> Chain 1: Iteration: 16 / 20 [ 80%]  (Sampling)
#> Chain 1: Iteration: 18 / 20 [ 90%]  (Sampling)
#> Chain 1: Iteration: 20 / 20 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.215 seconds (Warm-up)
#> Chain 1:                0.208 seconds (Sampling)
#> Chain 1:                0.423 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 2.12, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "compnet uses Stan under the hood. You may see warnings from Stan alongside, \nthis message. To deal with any warnings Stan might issue, \nPlease see the links provided in Stan's output, as well as the compnet website:\nhttps://kyle-rosenblad.github.io/compnet/"
fixedeff_violins(ex_compnet)

```
