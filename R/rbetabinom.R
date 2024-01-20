#' Simulate draws from a beta-binomial distribution
#'
#' @importFrom stats rbinom rbeta
#' @export
#' @param n Integer value: number of samples to draw.
#' @param size Integer value: number of trials for each sample.
#' @param prob Vector, matrix, or array of real values specifying the probability of success.
#' @param phi Real value: "over/under- dispersion parameter" for beta-binomial likelihood. In this
#'    parameterization, the latent beta variable has shape \eqn{\mu\phi} and scale
#'    \eqn{(1-\mu)\phi}, where \eqn{\mu} is the probability of success in the binomial distribution.
#' @return A vector of draws from the specified beta-binomial distribution.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- compnet(presabs=ex_presabs, spvars_dist_int=ex_traits, warmup=100, iter=200)
#' ex_compnet_pps <- postpredsamp(ex_compnet)
#'
rbetabinom <- function(n, size, prob, phi){
  samp <- stats::rbinom(n=n,
                        size=size,
                        prob=stats::rbeta(n=n,
                                          shape1=prob*phi,
                                          shape2=(1-prob)*phi))
  return(samp)
}
