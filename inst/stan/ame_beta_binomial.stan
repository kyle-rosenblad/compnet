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
  int K ; 
  int both[N] ; 
  int either[N] ; 
  real<lower=0> prior_intercept_scale ;
  real<lower=0> prior_betas_scale ;
  real<lower=0> prior_sigma_addeff_rate ;
  real<lower=0> prior_multi_cholesky_eta ;
  real<lower=0> prior_sigma_multi_shape ;
  real<lower=0> prior_sigma_multi_scale ;
  real<lower=0> prior_lambda_scale ;
  real<lower=0> prior_phi_rate ;
}
parameters{
  real<lower=0> phi ;

  real intercept ;
  
  vector[Xdy_cols] beta_dy ;
  vector[Xsp_cols] beta_sp ;
  
  real<lower=0> sigma_species_randint ; 
  vector[n_nodes] species_randint ; 
  
  cholesky_factor_corr[K] corr_multi_effects ; 
  vector<lower=0>[K] sigma_multi_effects ; 
  matrix[K, n_nodes] z_multi_effects ; 
  
  vector<lower=0>[K] lambda_diag ; 
}
transformed parameters{
    matrix[n_nodes, K] mean_multi_effects ; 
    mean_multi_effects = (diag_pre_multiply(
      sigma_multi_effects, corr_multi_effects) * z_multi_effects)' ;

	vector[N] xb ;
	xb = intercept + Xdy*beta_dy + XA*beta_sp + XB*beta_sp ;
	
	matrix[K, K] lambda ;
	lambda = diag_matrix(lambda_diag) ;
}
model{
  phi ~ exponential(prior_phi_rate);

  intercept ~ normal(0, prior_intercept_scale) ; 
  
  beta_dy ~ normal(0, prior_betas_scale);
  beta_sp ~ normal(0, prior_betas_scale);

  sigma_species_randint ~ exponential(prior_sigma_addeff_rate) ; 
  species_randint ~ normal(0, 1) ; 
  
  to_vector(z_multi_effects) ~ normal(0, 1) ; 
  corr_multi_effects ~ lkj_corr_cholesky(prior_multi_cholesky_eta) ; 
  sigma_multi_effects ~ gamma(prior_sigma_multi_shape, prior_sigma_multi_scale) ; 

  lambda_diag ~ normal(0, prior_lambda_scale) ; 

  vector[N] pboth ; 
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] +
		mean_multi_effects[spAid[n], 1:K] * lambda *
		(mean_multi_effects[spBid[n], 1:K]') ;
	
    pboth[n] = inv_logit(pboth[n]) ;
  }
    both ~ beta_binomial( either, pboth*phi, (1-pboth*phi) ) ;
}

generated quantities{
  vector[N] pboth ; 
  for(n in 1:N){
    pboth[n] = xb[n] +
		sigma_species_randint*species_randint[spAid[n]] +
		sigma_species_randint*species_randint[spBid[n]] +
		mean_multi_effects[spAid[n], 1:K] * lambda *
		(mean_multi_effects[spBid[n], 1:K]') ;
	
    pboth[n] = inv_logit(pboth[n]) ;
  }
}
