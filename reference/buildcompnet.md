# Network models of interspecific competition with presence-absence data

Network models of interspecific competition with presence-absence data

## Usage

``` r
buildcompnet(
  presabs,
  spvars_no_int = NULL,
  spvars_dist_int = NULL,
  spvars_multi_int = NULL,
  spvars_cat_no_int = NULL,
  spvars_cat_int = NULL,
  pairvars = NULL,
  rank = 0,
  family = "fnchypg",
  olre = TRUE,
  prior_intercept_scale = 5,
  prior_betas_scale = 5,
  prior_sigma_addeff_rate = 1,
  prior_sigma_olre_rate = 1,
  prior_lambda_scale = 5,
  warmup = 1000,
  iter = 2000,
  adapt_delta = 0.8
)
```

## Arguments

- presabs:

  Must be specified by user. Binary (0 or 1) presence-absence matrix
  with sites as rows and species as columns. Column names should be
  unique species names.

- spvars_no_int:

  A matrix or data frame, in which rows are species and columns are
  traits to be included in the model as additive species effects only,
  with no interaction term. Row names should be unique species names,
  and column names should be unique trait names.

- spvars_dist_int:

  A matrix or data frame, in which rows are species and columns are
  traits to be included in the model with a distance interaction term
  (i.e., the absolute value of the difference in trait values for each
  species pair). Row names should be unique species names, and column
  names should be unique trait names.

- spvars_multi_int:

  A matrix or data frame, in which rows are species and columns are
  traits to be included in the model with a multiplicative interaction
  term (i.e., the product of the trait values for each species pair).
  Row names should be unique species names, and column names should be
  unique trait names.

- spvars_cat_no_int:

  A matrix or data frame, in which rows are species and columns are
  categorical traits to be included in the model as additive species
  effects only, with no interaction term. Row names should be unique
  species names, and column names should be unique trait names.

- spvars_cat_int:

  A matrix or data frame, in which rows are species and columns are
  categorical traits to be included in the model with a binary "same or
  different" interaction term. Row names should be unique species names,
  and column names should be unique trait names.

- pairvars:

  A matrix or data frame, in which rows are species pairs and columns
  are pair-level traits that do not have a species-level analog, e.g.,
  phylogenetic distance. Column names should be unique trait names.
  There should also be two columns named "spAid" and "spBid" containing
  the unique names of the species in each pair.

- rank:

  Number of dimensions for the multiplicative latent factor term. Rank=0
  (the default) yields a model with no multiplicative term.

- family:

  Distribution family for regression model. Defaults to 'fnchypg', or
  Fisher's noncentral hypergeometric distribution. In this case, the
  quantity we are modeling, as a function of species traits, dyadic
  random effects, etc. is Mainali et al.'s (2022, Science Advances)
  "alpha" or "cooccurrence affinity" parameter in the fnchypg
  distribution. Link function is identity. 'binomial', also supported,
  is not as theoretically well justified– see Mainali et al. (2022,
  Science Advances)–but, in simulations, produces qualitatively similar
  results, runs faster, and more consistently avoids problems like
  overdispersion. In this case, we are modeling p, the probability that
  both species co-occur at a given site, given that at least one is
  present. Link function is logit. May be useful for pilot analyses.

- olre:

  Logical variable indicating whether to include observation-level
  random intercepts drawn from a univariate normal distribution. Only
  relevant when family = 'fnchypg'. Defaults to TRUE, as these models
  can often exhibit fit problems which may be mitigated by the random
  effects.

- prior_intercept_scale:

  Scale parameter for mean-zero Gaussian prior on the intercept term for
  the linear predictor.

- prior_betas_scale:

  Scale parameter for mean-zero Gaussian priors on the coefficients of
  fixed effect terms in the linear predictor.

- prior_sigma_addeff_rate:

  Rate parameter for exponential prior on the scale of the species-level
  Gaussian random effects (i.e., "row and column effects").

- prior_sigma_olre_rate:

  Rate parameter for exponential prior on the scale of the dyad-level
  Gaussian random effects.

- prior_lambda_scale:

  Scale parameter for mean-zero Gaussian prior on diagonal values of
  Lambda, the matrix that determines how different species' values of
  the latent factors interact in the linear predictor.

- warmup:

  Number of warmup iterations for Stan.

- iter:

  Number of posterior sampling iterations for Stan.

- adapt_delta:

  A parameter that tunes Stan's posterior sampling algorithm. Increasing
  closer to 1 can help avoid divergent transitions.

## Value

Object of class "compnet", which is a list containing the stanfit model
object, a named list of posterior samples for all model parameters, a
data frame containing all input variables for the model, a matrix of
dyadic X variables, a matrix of X variables pertaining to species A in
each pair, a matrix of X variables pertaining to species B in each pair,
and–when relevant–a matrix of means and standard deviations for the
input trait data before centering and scaling.

## Details

This function uses Stan, as accessed through the rstan package, to build
a network regression model of interspecific competitive niche
differentiation in a Bayesian framework. This function is designed to
test the hypothesis that species are more likely to co-occur with other
species that are functionally dissimilar. Functional dissimilarity can
be represented directly by traits or by a proxy like phylogenetic
distance. The core input data are a species-by-site presence-absence
matrix and one or more species-level traits (e.g., plant leaf size) or
pair-level traits (e.g., phylogenetic distance). Units of analysis are
species pairs. The response variable follows a binomial distribution.
The number of trials is the number of sites occupied by at least one
species in a pair, and the number of successes is the number of sites
occupied by both species. If species-level traits are used, each trait
can be non-interacting (i.e., there is no interaction term between
species A's trait value and species B's), interacting via a typical
multiplicative term, or interacting via an absolute value difference
(i.e., "distance") term. The interaction terms are key to the core
hypothesis. If competitive niche differentiation is occurring, then the
probability of co-occurrence is expected to increase with trait or
phylogenetic distance between species A and B, or with the product of
their trait values, if a multiplicative interaction is specified instead
of a distance interaction. Random effects are used to account for
additive species-level dependencies and, optionally, higher-order
dependencies involving multiple species (e.g., "the enemy of my enemy is
my friend"). Higher-order dependencies are modeled using a number of
latent variable specified by "rank". For more details on the random
effects, see Hoff, P. (2021) Additive and multiplicative effects network
models. Stat. Sci. 36, 34–50.

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
#> Chain 1: Gradient evaluation took 0.000399 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.99 seconds.
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
#> Chain 1:  Elapsed Time: 0.369 seconds (Warm-up)
#> Chain 1:                0.211 seconds (Sampling)
#> Chain 1:                0.58 seconds (Total)
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
```
