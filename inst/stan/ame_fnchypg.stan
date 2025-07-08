functions {
  real fnchypg_lpmf(int y, int mA, int mB, int n, real alpha) {
    // set up possible y values
    int L = max(0, mA + mB - n);
    int U = min(mA, mB);

    // set up useful objects for looping through pmf calculations
    int len = U - L + 1;
    array[len] int x = linspaced_int_array(len, L, U);
    vector[len] log_probs;

    // compute lpm
    for (i in 1:len) {
      log_probs[i] = lchoose(mA, x[i]) + lchoose(n - mA, mB - x[i]) + x[i] * log(alpha);
    }

    return log_probs[y - L + 1] - log_sum_exp(log_probs);
  }
}

data{
  int n_nodes; // number of species
  int N; // number of dyads
  int<lower=1> n_sites; // number of sites
  array[N] int spAid; // ID number for species A in each dyad
  array[N] int spBid; // ID number for species B in each dyad
  int Xdy_cols; // number of columns in the matrix of dyad-level predictor variables
  int Xsp_cols;  // number of columns in the matrix of species-level predictor variables
  matrix[N,Xdy_cols] Xdy; // matrix of dyad-level predictor variables
  matrix[N,Xsp_cols] XA;  // matrix of species-level predictor variables with values for species A
  matrix[N,Xsp_cols] XB;  // matrix of species-level predictor variables with values for species A
  int K; // number of dimensions for the latent factor model
  array[N] int both; // for each dyad, the number of sites occupied by both species in each dyad
  array[n_nodes] int sp_occ;  // total site occupancy for each species
  real<lower=0> prior_intercept_scale; // SD for normal prior on the intercept
  real<lower=0> prior_betas_scale; // SD for normal prior on the coefficients of predictor variables
  real<lower=0> prior_sigma_addeff_rate; // rate for exponential prior on the sd of the species-level random effects (i.e., "additive" effects)
  real<lower=0> prior_lambda_scale; // SD for normal prior on diagonal values of the lambda matrix, which scales the strength of the effect of each latent dimension on cooccurrence probabilities
}

parameters{
  real intercept; // intercept for linear predictor
  vector[Xdy_cols] beta_dy; // coefficients for dyad-level terms in linear predictor
  vector[Xsp_cols] beta_sp; // coefficients for species-level terms in linear predictor
  real<lower=0> sigma; // sd of species-level random effects (i.e., "additive" effects)
  vector[n_nodes] a_raw; // species-level random effects (i.e., "additive" effects)
  matrix[n_nodes, K] U_raw; // multiplicative effect values for each species across each latent dimension in latent factor model
  ordered[K] lambda_diag; // diagnoal values of the lambda matrix, which controls the strength and direction of multiplicative effects
}

transformed parameters{
    // scale species-level random intercepts from non-centered parameterization
    vector[n_nodes] a = sigma * a_raw;

  // Set the first element in each latent factor as positive to avoid sampling
  // issues with reflectional invariance, which we can think of pragmatically
  // as a form of nonidentifiability. i.e., the multiplicative inverse ("times -1") of any given set of
  // latent factor values produces an identical likelihood, and we don't want the sampler
  // to waste time toggling between functionally equivalent mirror-image configurations of latent factors.
  matrix[n_nodes, K] U;
  for (k in 1:K) {
    U[1, k] = fabs(U_raw[1, k]);  // get absolute value of first node
    for (i in 2:n_nodes) {
      U[i, k] = U_raw[i, k];      // leave other nodes alone
    }
  }
}

model{
  // set priors
  intercept ~ normal(0, prior_intercept_scale);
  beta_dy ~ normal(0, prior_betas_scale);
  beta_sp ~ normal(0, prior_betas_scale);

  // Use non-centered parameterization for species random effects.
  // This is a trick to make the sampler work better.
  // Instead of initially defining species_randint ~ normal(0, sigma_species_randint),
  // we z-score the random effects and then multiply by sigma later.
  a_raw ~ normal(0, 1);
  sigma ~ exponential(prior_sigma_addeff_rate);

  // Raw latent factors are all standard normal, as shown below, and then
  // the diagonal values of lambda determine the strength and direction
  // of the multiplicative interaction between species for each latent
  // dimension.
  for (i in 1:n_nodes) {
    U_raw[i] ~ normal(0, 1);
  }
  // prior for the diagonals of lambda, as described above
  lambda_diag ~ normal(0, prior_lambda_scale);

  // likelihood
  for (i in 1:N) {
    real eta = intercept +
               dot_product(Xdy[i], beta_dy) +
               dot_product(XA[i], beta_sp) +
               dot_product(XB[i], beta_sp) +
               a[spAid[i]] + a[spBid[i]] +
               U[spAid[i], 1:K] * diag_matrix(lambda_diag) *
        (U[spBid[i], 1:K]');
    real alpha = exp(eta); // alpha or "affinity" for Fisher's noncentral hypergeometric distribution
    int mA_i = sp_occ[spAid[i]]; // get occurrence frequency for species A in dyad i
    int mB_i = sp_occ[spBid[i]]; // get occurrence frequency for species B in dyad i
    target += fnchypg_lpmf(both[i] | mA_i, mB_i, n_sites, alpha);
  }
}
generated quantities{
  vector[N] eta;
  for(i in 1:N){
    eta[i] = intercept +
               dot_product(Xdy[i], beta_dy) +
               dot_product(XA[i], beta_sp) +
               dot_product(XB[i], beta_sp) +
               a[spAid[i]] + a[spBid[i]] +
               U[spAid[i], 1:K] * diag_matrix(lambda_diag) *
        (U[spBid[i], 1:K]');
  }
  vector[N] latfacterm;
  for(i in 1:N){
    latfacterm[i] = U[spAid[i], 1:K] * diag_matrix(lambda_diag) *
        (U[spBid[i], 1:K]');
  }
}
