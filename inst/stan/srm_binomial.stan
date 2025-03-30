data{
  int n_nodes; // number of species
  int N; // number of dyads
  array[N] int spAid;  // ID number for species A in each dyad
  array[N] int spBid;  // ID number for species B in each dyad
  int Xdy_cols; // number of columns in the matrix of dyad-level predictor variables
  int Xsp_cols;  // number of columns in the matrix of species-level predictor variables
  matrix[N,Xdy_cols] Xdy; // matrix of dyad-level predictor variables
  matrix[N,Xsp_cols] XA;  // matrix of species-level predictor variables with values for species A
  matrix[N,Xsp_cols] XB;  // matrix of species-level predictor variables with values for species A
  array[N] int both; // for each dyad, the number of sites occupied by both species in each dyad
  array[N] int either; // for each dyad, the number of sites occupied by at least one of the two species
  real<lower=0> prior_intercept_scale; // SD for normal prior on the intercept
  real<lower=0> prior_betas_scale; // SD for normal prior on the coefficients of predictor variables
  real<lower=0> prior_sigma_addeff_rate; // rate for exponential prior on the sd of the species-level random effects (i.e., "additive" effects)
}
parameters{
  real intercept; // intercept for linear predictor

  vector[Xdy_cols] beta_dy; // coefficients for dyad-level terms in linear predictor
  vector[Xsp_cols] beta_sp; // coefficients for species-level terms in linear predictor

  real<lower=0> sigma_species_randint; // sd of species-level random effects (i.e., "additive" effects)
  vector[n_nodes] species_randint; // species-level random effects (i.e., "additive" effects)
}
transformed parameters{
  // bundle all fixed effects in "xb" for convenience
  vector[N] xb;
  xb = intercept + Xdy*beta_dy + XA*beta_sp + XB*beta_sp;
}
model{
  // set priors
  intercept ~ normal(0, prior_intercept_scale);
  beta_dy ~ normal(0, prior_betas_scale);
  beta_sp ~ normal(0, prior_betas_scale);
  sigma_species_randint ~ exponential(prior_sigma_addeff_rate);
  // Use non-centered parameterization for species random effects.
  // This is a trick to make the sampler work better.
  // Instead of initially defining species_randint ~ normal(0, sigma_species_randint),
  // we z-score the random effects and then multiply by sigma later.
  species_randint ~ normal(0, 1);

  // put it all together with a logit link function to the binomial likelihood:
  vector[N] pboth ;
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] ;

    pboth[n] = inv_logit(pboth[n]) ;
  }
    both ~ binomial( either , pboth ) ;
}
// Get some convenient quantities for model diagnostics:
generated quantities{
  vector[N] pboth ;
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] ;

    pboth[n] = inv_logit(pboth[n]) ;
  }
}
