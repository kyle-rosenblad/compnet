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
#' bbdraws <- rbetabinom(n=100, size=10, prob=0.3, phi=2)
#'
rbetabinom <- function(n, size, prob, phi){
  samp <- stats::rbinom(n=n,
                        size=size,
                        prob=stats::rbeta(n=n,
                                          shape1=prob*phi,
                                          shape2=(1-prob)*phi))
  return(samp)
}
