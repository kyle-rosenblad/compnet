The goal: disentangle the roles of multiple, possibly correlated species
traits and phylogenetic relationships in shaping species co-occurrence.
‘compnet’ uses Bayesian dyadic regression models to quantify these
effects with presence-absence data and species-level trait data (e.g.,
plant leaf size) and/or pair-level trait data (e.g., phylogenetic
distance). Specialized random effects and latent variables are used to
deal with the types of non-independence found in network data. First,
follow the rstan installation instructions:
(<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>). Then,
in R, run install.packages(“devtools”), then
devtools::install_github(“kyle-rosenblad/compnet”).
