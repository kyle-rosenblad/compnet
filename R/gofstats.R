#' Report quantiles of gofstats values for observed data relative to posterior predictive distribution
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param thin Logical value indicating whether to take a random subsample of posterior draws.
#' @param thin_to Logical value indicating the size of the subsample to take if thin=TRUE.
#' @return A named vector containing quantiles on the interval \eqn{[0,1]} for: 1- the standard deviation
#'    of row means, and 2- The triadic dependency metric used by Hoff, Fosdick, & Volfovsky's "amen"
#'    package. These values represent the proportion of the posterior predictive simulation that were
#'    less than the value for the observed data.
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
#' ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#'
#' gofstats(ex_compnet)
#'

gofstats <- function(mod,
                     thin=T,
                     thin_to=300){
  d <- mod$d
  N_species <- length(unique(c(d$spAid, d$spBid)))

  if(mod$family=='fnchypg'){

    dsquare <- matrix(NA, nrow = N_species, ncol = N_species)
    for(i in 1:nrow(dsquare)){
      tempvec <- subset(d, d$spAid==i | d$spBid==i)$both
      tempvec <- append(tempvec, NA, after = i - 1)
      dsquare[i, ] <- tempvec
    }

    gofstats_obs <- compute_gofstats(dsquare)
    gofstats_sim <- data.frame(sd.rowmean = NA, cycle.dep = NA)

    # Get posterior predictive matrix from fnchypg model
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
        tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$ppred
        tempvec <- append(tempvec, NA, after = i - 1)
        dsquaretemp[i,] <- tempvec
      }
      gofstats_sim[j,] <- compute_gofstats(dsquaretemp)
    }

  }



  if(mod$family=='binomial'){
    dsquare <- matrix(NA, nrow = N_species, ncol = N_species)
    for(i in 1:nrow(dsquare)){
      tempvec <- subset(d, d$spAid==i | d$spBid==i)$both
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
        tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$ppred
        tempvec <- append(tempvec, NA, after = i - 1)
        dsquaretemp[i,] <- tempvec
      }
      gofstats_sim[j,] <- compute_gofstats(dsquaretemp)
    }
  }


  p.sd.rowmeans <- mean(gofstats_obs[1] > gofstats_sim[, 1])
  p.cycle.dep   <- mean(gofstats_obs[2] > gofstats_sim[, 2])
  out <- c(p.sd.rowmeans, p.cycle.dep)
  names(out) <- c("p.sd.rowmeans", "p.cycle.dep")
  return(out)
}
