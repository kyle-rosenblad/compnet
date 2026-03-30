# Make an automated plot of the interactive effect of two species' values of the same trait.

Make an automated plot of the interactive effect of two species' values
of the same trait.

## Usage

``` r
plot_interaction(
  mod,
  xvar,
  xlabel,
  orig.scale = TRUE,
  intlevels = c(0.05, 0.5, 0.95),
  ci_width = 0.95,
  grid_size = 100,
  thin = TRUE,
  thin_to = 100
)
```

## Arguments

- mod:

  An object of class "compnet" created by the buildcompnet() function.

- xvar:

  Character string for the name of the trait to be used. Must match the
  trait name in the input data used to build the model.

- xlabel:

  Optional character string to replace xvar when plotting.

- orig.scale:

  Logical value indicating whether to back-transform trait data to the
  original scale (TRUE) or leave them with mean zero and unit variance
  (FALSE).

- intlevels:

  Vector of real values on the interval \\\[0,1\]\\ indicating what
  levels of the x variable to condition on for species B when plotting
  species A's mean response.

- ci_width:

  A real number (0,1) describing the desired widths of credible bands.
  Defaults to 0.95.

- grid_size:

  A positive integer defining the number of discrete steps to use in
  approximating the shape of mean prediction curves and credible bands.
  Defaults to 100.

- thin:

  Logical value determining whether to use a random subsample of the
  full posterior sample.

- thin_to:

  Integer value determining how many random samples to draw from the
  full posterior sample if thin=TRUE.

## Value

A ggplot2 graphic.

## Examples

``` r
data(ex_presabs)
data(ex_traits)

# Quick demo run. Will prompt warnings.
# Run with default warmup and iter for good posterior sampling.
ex_compnet <- buildcompnet(presabs=ex_presabs,
spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20, family='binomial')
#> 
#> SAMPLING FOR MODEL 'srm_binomial' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000114 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.14 seconds.
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
#> Chain 1:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 1:                0.05 seconds (Sampling)
#> Chain 1:                0.078 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
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
plotdata <- plot_interaction(ex_compnet, xvar="ndtrait")
```
