#' Simulate posterior predictive samples from a compnet model
#'
#' @importFrom BiasedUrn rFNCHypergeo
#' @importFrom stats rbinom
#' @export
#' @param mod Object of class "compnet", which is created by the buildcompnet() function.
#' @return A matrix of posterior predictive samples with a row for each observation and a column for each sample.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- buildcompnet(presabs=ex_presabs,
#' spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#' ex_compnet_pps <- postpredsamp(ex_compnet)
#'

postpredsamp <- function(mod){

  if(mod$family=='fnchypg'){
    d <- mod$d
    sp_occ <- mod$datalist$sp_occ
    spAocc <- data.frame(spAid_orig=names(sp_occ),
                         freqA=sp_occ)
    spBocc <- data.frame(spBid_orig=names(sp_occ),
                         freqB=sp_occ)
    d <- merge(d, spAocc, by="spAid_orig", all.x=TRUE, all.y=FALSE, sort=FALSE)
    d <- merge(d, spBocc, by="spBid_orig", all.x=TRUE, all.y=FALSE, sort=FALSE)

    postpredsamp <- t(mod$stanmod_samp$alpha) # make this a draw from the distribution rather than just alpha

    for(i in 1:nrow(postpredsamp)){
      postpredsamp[i,] <- sapply(exp(postpredsamp[i,]), BiasedUrn::rFNCHypergeo, nran=1,
                                 m1=d$freqA[i], m2=mod$datalist$n_sites-d$freqA[i], n=d$freqB[i])

    }

  }


  if(mod$family=='binomial'){
    postpredsamp <- t(mod$stanmod_samp$pboth)
    for(i in 1:nrow(postpredsamp)){
      postpredsamp[i,] <- stats::rbinom(n=ncol(postpredsamp),
                                        size=mod$d[i, "either"],
                                        prob=postpredsamp[i,])
    }

  }


  return(postpredsamp)
}
