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
  }
    both ~ binomial( either , pboth ) ;
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
