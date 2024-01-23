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
#' ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits, warmup=10, iter=20)
#'
#' gofstats(ex_compnet)
#'

gofstats <- function(mod,
                     thin=T,
                     thin_to=300){
  d <- mod$d
  d$pboth <- d$both/d$either
  dsquare <- matrix(NA, nrow=length(unique(c(d$spAid, d$spBid))), ncol=length(unique(c(d$spAid, d$spBid))))
  for(i in 1:nrow(dsquare)){
    tempvec <- subset(d, d$spAid==i | d$spBid==i)$pboth
    if(i==1){
      tempvec <- c(NA, tempvec)
    }
    if(i>1 & i<length(unique(c(d$spAid, d$spBid)))){
      tempvec <- c(tempvec[c(1:(i-1))], NA, tempvec[c(i:length(tempvec))])
    }
    if(i==length(unique(c(d$spAid, d$spBid)))){
      tempvec <- c(tempvec, NA)
    }
    dsquare[i,] <- tempvec
  }

  gofstats_obs <- compute_gofstats(dsquare)
  gofstats_sim <- data.frame(sd.rowmean=NA,
                             cycle.dep=NA)

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
    dtemp$pboth2 <- dtemp$ppred/dtemp$either
    dsquaretemp <- matrix(NA, nrow=length(unique(c(d$spAid, d$spBid))), ncol=length(unique(c(d$spAid, d$spBid))))

    for(i in 1:nrow(dsquaretemp)){
      tempvec <- subset(dtemp, dtemp$spAid==i | dtemp$spBid==i)$pboth2
      if(i==1){
        tempvec <- c(NA, tempvec)
      }
      if(i>1 & i<length(unique(c(d$spAid, d$spBid)))){
        tempvec <- c(tempvec[c(1:(i-1))], NA, tempvec[c(i:length(tempvec))])
      }
      if(i==length(unique(c(d$spAid, d$spBid)))){
        tempvec <- c(tempvec, NA)
      }
      dsquaretemp[i,] <- tempvec
    }
    gofstats_sim[j,] <- compute_gofstats(dsquaretemp)
  }
  p.sd.rowmeans <- sum(as.numeric(gofstats_obs[1]>gofstats_sim[,1]))/nrow(gofstats_sim)
  p.cycle.dep <- sum(as.numeric(gofstats_obs[2]>gofstats_sim[,2]))/nrow(gofstats_sim)
  out <- c(p.sd.rowmeans, p.cycle.dep)
  names(out) <- c("p.sd.rowmeans", "p.cycle.dep")
  return(out)
}
