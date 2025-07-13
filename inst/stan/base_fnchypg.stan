functions {
  real fnchypg_lpmf(int both, int Apres, int Aabs, int Bpres, real alpha) {
    // Set up the region of supported values for Fisher's noncentral hypergeometric distribution,
    // given the number of sites occupied by each species and the total number of sites sampled.
    // For example, if each species only occupies 1 site, there can't be more than 1 site where they co-occur.
    // Similarly, if each species occupies 99 sites, they must co-occur at at least 98 sites.
    int lowerbound = max(0, Bpres - Aabs);
    int upperbound = min(Apres, Bpres);
    int supportrange = upperbound - lowerbound + 1;

    // Compute log(probability mass). This is essentially an exercise in combinatorics,
    // whereby we explore all the ways to get a given number of sites with co-occurrence,
    // given how many sites are occupied by each species, and how many total sampling sites
    // there are. If you are familiar with "biased urn" problems, it can be helpful to think
    // of it that way. For each species pair, the total number of balls in the urn represents
    // the total number of sampling sites. The red balls are sites where species A is present,
    // and the white balls are sites where species A is absent. The number of balls drawn from
    // the urn represents the number of sites where species B is present. The alpha parameter
    // introduces bias to the urn, based on species' propensity to cooccur with each other
    // more or less often than they would in an unbiased urn scenario.
    vector[supportrange] log_probs;
    for (i in 1:supportrange) {
      int x = lowerbound + i - 1;
      real log_choose_term1 = lchoose(Apres, x);
      real log_choose_term2 = lchoose(Aabs, Bpres - x);
      log_probs[i] = log_choose_term1 + log_choose_term2 + x * alpha;
    }
    // Return log probability of observed value, normalized to the full probability mass function
    int index = both - lowerbound + 1;
    return log_probs[index] - log_sum_exp(log_probs);
  }
}
data{
  int n_nodes; // number of species
  int N; // number of dyads
  int<lower=1> n_sites; // number of sites
  array[N] int spAid; // ID number for species A in each dyad
  array[N] int spBid; // ID number for species B in each dyad
  array[N] int both; // for each dyad, the number of sites occupied by both species in each dyad
  array[n_nodes] int sp_occ;  // total site occupancy for each species
}

parameters {
  real mu; // mean of alpha distribution
  real<lower=0> sigma; // sd of alpha distribution
  vector[N] alpha_raw; // standardized alpha values
}

transformed parameters {
  vector[N] alpha = mu + sigma * alpha_raw; // non-centered parameterization
}
model {
  // Priors
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  alpha_raw ~ normal(0,1); // implies alpha ~ normal(mu, sigma)

  // Likelihood
  for (i in 1:N) {
    int spApres_i = sp_occ[spAid[i]]; // number of sites occupied by species A in dyad i
    int spAabs_i = n_sites - sp_occ[spAid[i]]; // number of sites NOT occupied by species A in dyad i
    int spBpres_i = sp_occ[spBid[i]]; // number of sites occupied by species B in dyad i

    target += fnchypg_lpmf(both[i] | spApres_i, spAabs_i, spBpres_i, alpha[i]);
  }
}
