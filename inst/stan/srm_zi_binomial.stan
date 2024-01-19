functions {
  real zero_inflated_binomial_lpmf(int y, int trials,
                                   real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  real zero_inflated_binomial_logit_lpmf(int y, int trials,
                                         real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  real zero_inflated_binomial_blogit_lpmf(int y, int trials,
                                          real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  real zero_inflated_binomial_blogit_logit_lpmf(int y, int trials,
                                                real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  real zero_inflated_binomial_lccdf(int y, int trials, real theta, real zi) {
    return bernoulli_lpmf(0 | zi) + binomial_lccdf(y | trials, theta);
  }
  real zero_inflated_binomial_lcdf(int y, int trials, real theta, real zi) {
    return log1m_exp(zero_inflated_binomial_lccdf(y | trials, theta, zi));
  }
}
data{
  int n_nodes ;
  int N ; 
  int spAid[N] ;  
  int spBid[N] ;  
  int Xdy_cols ; 
  int Xsp_cols ;  
  matrix[N,Xdy_cols] Xdy ; 
  matrix[N,Xsp_cols] XA ;  
  matrix[N,Xsp_cols] XB ;  
  int both[N] ; 
  int either[N] ; 
  real<lower=0> prior_intercept_scale ;
  real<lower=0> prior_betas_scale ;
  real<lower=0> prior_sigma_addeff_rate ;
}
parameters{
  real<lower=0, upper=1> zi;

  real intercept ;

  vector[Xdy_cols] beta_dy ;
  vector[Xsp_cols] beta_sp ;
  
  real<lower=0> sigma_species_randint ; 
  vector[n_nodes] species_randint ;   
}
transformed parameters{
	vector[N] xb ;
	xb = intercept + Xdy*beta_dy + XA*beta_sp + XB*beta_sp ;
}
model{
  intercept ~ normal(0, prior_intercept_scale) ; 
  
  beta_dy ~ normal(0, prior_betas_scale);
  beta_sp ~ normal(0, prior_betas_scale);

  sigma_species_randint ~ exponential(prior_sigma_addeff_rate) ; 
  species_randint ~ normal(0, 1) ; 
  
  vector[N] pboth ; 
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] ;
	
    pboth[n] = inv_logit(pboth[n]) ;
	
	target += zero_inflated_binomial_blogit_lpmf(both[n] | either[n], pboth[n], zi) ;

  }
}

generated quantities{
  vector[N] pboth ; 
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] ;
	
    pboth[n] = inv_logit(pboth[n]) ;
  }
}
