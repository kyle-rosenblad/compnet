# Summarizing compnet model output

Summarizing compnet model output

## Usage

``` r
summarize_compnet(mod, ci_width = 0.95)
```

## Arguments

- mod:

  Object of class "compnet", which is created by the buildcompnet()
  function.

- ci_width:

  A real number (0,1) of the desired interval width. Defaults to 0.95.

## Value

A data frame summarizing means and credible intervals for standardized
effect sizes of fixed effects.

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
#> Chain 1: Gradient evaluation took 0.000407 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 4.07 seconds.
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
#> Chain 1:  Elapsed Time: 0.389 seconds (Warm-up)
#> Chain 1:                0.256 seconds (Sampling)
#> Chain 1:                0.645 seconds (Total)
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

ex_compnet_summ <- summarize_compnet(ex_compnet)
ex_compnet_summ
#>                    Mean       2.5%      97.5%
#> intercept    -2.0488194 -2.8486059 -1.2478814
#> ndtrait_dist  0.9438650  0.7082897  1.2155203
#> ndtrait_sp    0.2341691 -0.3980288  0.6835962
```
