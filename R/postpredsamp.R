#' Simulate posterior predictive samples from a compnet model
#'
#' @importFrom stats rbinom
#' @importFrom BiasedUrn rFNCHypergeo
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
#' ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#' ex_compnet_pps <- postpredsamp(ex_compnet)
#'

postpredsamp <- function(mod){

  if(mod$family=='fnchypgm'){
    postpredsamp <- t(mod$stanmod_samp$alpha)
    #for(i in 1:nrow(postpredsamp)){
    #  postpredsamp[i,] <- BiasedUrn::rFNCHypergeo(nran=1,
    #                                              m1=mod$d$datalist$
    #                                              odds=exp(postpredsamp[i,]))
    #}

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
