#' Report quantiles of gofstats values for observed data relative to posterior predictive distribution
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param thin Logical value indicating whether to take a random subsample of posterior draws.
#' @param thin_to Logical value indicating the size of the subsample to take if thin=TRUE.
#' @return A named vector containing quantiles on the interval \eqn{[0,1]} for: 1- the standard deviation
#'    of row means, and 2- The triadic dependency metric used by Hoff, Fosdick, & Volfovsky's "amen"
#'    package. These values represent the proportion of the posterior predictive simulation that were
#'    less than the value for the reference comparison. For models fit with Fisher's noncentral hypergeometric.
#'    distribution, this reference comparison is the distribution of alpha from a model fit with vague priors
#'    and no random effects, fixed effects, etc., to represent the degree of dyadic and triadic nonindependence
#'    in alpha in the original data as faithfully as possible. (Using posterior predictive draws of "both" does
#'    not work well for technical reasons having to do with the complex nature of the FNCHYPG likelihood). For
#'    models fit with the binomial distribution, the reference comparison is the ratio of "both" to "either" in the original data.
#' @details This function can be used to assess whether species-level and higher-order dependencies
#'    in the data are represented adequately by the model structure. Extreme output values indicate
#'    there may be a problem. In these cases, it may help to use different fixed effect predictors,
#'    or to increase the model rank.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- buildcompnet(presabs=ex_presabs,
#' spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#'
#' gofstats(ex_compnet)
#'

gofstats <- function(mod,
                     thin=T,
                     thin_to=300){
  d <- mod$d
  N_species <- length(unique(c(d$spAid, d$spBid)))

  if(mod$family=='fnchypg'){
    datalist=mod$datalist[c("n_nodes",
                             "N",
                             "sp_occ",
                             "spAid",
                             "spBid",
                             "both",
                             "n_sites")]
    cat("Fitting base model for comparison with full model")
    basemod <- rstan::sampling(stanmodels$base_fnchypg,
                               data=datalist,
                               cores=2,
                               chains=2,
                               warmup=1000,
                               iter=2000,
                               verbose=F,
                               control=list(adapt_delta=0.8))

    gofstats_base <- data.frame(sd.rowmean = NA, cycle.dep = NA)
    base_alpha <- t(rstan::extract(basemod)$alpha)
    if(thin==T){
      base_alpha <- base_alpha[, sample(1:ncol(base_alpha), size=min(thin_to, ncol(base_alpha)), replace=FALSE)]
    }
    for(j in 1:ncol(base_alpha)){
      dtemp <- d
      dtemp$base_alpha <- base_alpha[,j]
      dsquaretemp <- matrix(NA, nrow=N_species, ncol=N_species)

      for(i in 1:nrow(dsquaretemp)){
        tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$base_alpha
        tempvec <- append(tempvec, NA, after = i - 1)
        dsquaretemp[i,] <- tempvec
      }
      gofstats_base[j,] <- compute_gofstats(dsquaretemp)
    }

    gofstats_full <- data.frame(sd.rowmean = NA, cycle.dep = NA)
    full_alpha <- t(mod$stanmod_samp$alpha)
    if(thin==T){
      full_alpha <- full_alpha[, sample(1:ncol(full_alpha), size=min(thin_to, ncol(full_alpha)), replace=FALSE)]
    }
    for(j in 1:ncol(full_alpha)){
      dtemp <- d
      dtemp$full_alpha <- full_alpha[,j]
      dsquaretemp <- matrix(NA, nrow=N_species, ncol=N_species)

      for(i in 1:nrow(dsquaretemp)){
        tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$full_alpha
        tempvec <- append(tempvec, NA, after = i - 1)
        dsquaretemp[i,] <- tempvec
      }
      gofstats_full[j,] <- compute_gofstats(dsquaretemp)
    }
    p.sd.rowmeans <- mean(outer(gofstats_base$sd.rowmean, gofstats_full$sd.rowmean, ">"))
    p.cycle.dep   <- mean(outer(gofstats_base$cycle.dep, gofstats_full$cycle.dep, ">"))
    out <- c(p.sd.rowmeans, p.cycle.dep)
    names(out) <- c("p.sd.rowmeans", "p.cycle.dep")
  }



  if(mod$family=='binomial'){
    dsquare <- matrix(NA, nrow = N_species, ncol = N_species)
    for(i in 1:nrow(dsquare)){
      tempvec <- subset(d, d$spAid==i | d$spBid==i)$both / subset(d, d$spAid==i | d$spBid==i)$either
      tempvec <- append(tempvec, NA, after = i - 1)
      dsquare[i, ] <- tempvec
    }

    gofstats_obs <- compute_gofstats(dsquare)
    gofstats_sim <- data.frame(sd.rowmean = NA, cycle.dep = NA)

    ppred <- postpredsamp(mod)
    if(thin==T){
      ppred <- ppred[, sample(1:ncol(ppred), size=min(thin_to, ncol(ppred)), replace=FALSE)]
    }
    odens <- ncol(ppred)/4
    cat("Approx. completion", fill=T)
    for(j in 1:ncol(ppred)){
      if(j%%odens==0){
        cat(25*j/odens, "%", sep="", fill=T)
      }
      dtemp <- d
      dtemp$ppred <- ppred[,j]
      dsquaretemp <- matrix(NA, nrow=N_species, ncol=N_species)

      for(i in 1:nrow(dsquaretemp)){
        tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$ppred / subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$either
        tempvec <- append(tempvec, NA, after = i - 1)
        dsquaretemp[i,] <- tempvec
      }
      gofstats_sim[j,] <- compute_gofstats(dsquaretemp)
    }
    p.sd.rowmeans <- mean(gofstats_obs[1] > gofstats_sim[, 1])
    p.cycle.dep   <- mean(gofstats_obs[2] > gofstats_sim[, 2])
    out <- c(p.sd.rowmeans, p.cycle.dep)
    names(out) <- c("p.sd.rowmeans", "p.cycle.dep")
  }

  return(out)
}
