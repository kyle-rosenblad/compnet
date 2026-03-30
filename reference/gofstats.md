# Report quantiles of gofstats values for observed data relative to posterior predictive distribution

Report quantiles of gofstats values for observed data relative to
posterior predictive distribution

## Usage

``` r
gofstats(mod, thin = T, thin_to = 300)
```

## Arguments

- mod:

  An object of class "compnet" created by the buildcompnet() function.

- thin:

  Logical value indicating whether to take a random subsample of
  posterior draws.

- thin_to:

  Logical value indicating the size of the subsample to take if
  thin=TRUE.

## Value

A named vector containing quantiles on the interval \\\[0,1\]\\ for: 1-
the standard deviation of row means, and 2- The triadic dependency
metric used by Hoff, Fosdick, & Volfovsky's "amen" package. These values
represent the proportion of the posterior predictive simulation that
were less than the value for the reference comparison. For models fit
with Fisher's noncentral hypergeometric. distribution, this reference
comparison is the distribution of alpha from a model fit with vague
priors and no random effects, fixed effects, etc., to represent the
degree of dyadic and triadic nonindependence in alpha in the original
data as faithfully as possible. (Using posterior predictive draws of
"both" does not work well for technical reasons having to do with the
complex nature of the FNCHYPG likelihood). For models fit with the
binomial distribution, the reference comparison is the ratio of "both"
to "either" in the original data.

## Details

This function can be used to assess whether species-level and
higher-order dependencies in the data are represented adequately by the
model structure. Extreme output values indicate there may be a problem.
In these cases, it may help to use different fixed effect predictors, or
to increase the model rank.

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
#> Chain 1: Gradient evaluation took 0.000409 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 4.09 seconds.
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
#> Chain 1:  Elapsed Time: 0.214 seconds (Warm-up)
#> Chain 1:                0.212 seconds (Sampling)
#> Chain 1:                0.426 seconds (Total)
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

gofstats(ex_compnet)
#> Fitting base model for comparison with full model
#> SAMPLING FOR MODEL 'base_fnchypg' NOW (CHAIN 1).
#> 
#> SAMPLING FOR MODEL 'base_fnchypg' NOW (CHAIN 2).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000568 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 5.68 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000698 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 6.98 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 7.01 seconds (Warm-up)
#> Chain 1:                4.557 seconds (Sampling)
#> Chain 1:                11.567 seconds (Total)
#> Chain 1: 
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 7.007 seconds (Warm-up)
#> Chain 2:                4.562 seconds (Sampling)
#> Chain 2:                11.569 seconds (Total)
#> Chain 2: 
#> p.sd.rowmeans   p.cycle.dep 
#>     0.4560000     0.9683333 
```
