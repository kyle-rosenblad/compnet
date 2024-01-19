#' Inverse logit function
#'
#' @export
#' @param x A numeric value, or a vector or array of numeric values.
#' @return The inverse-logit(s) of the supplied numeric value(s).
#' @examples
#'
#' expit(-1)
#'
#' exvec <- c(-2:2)
#' expit(exvec)
#'
#' exmat <- matrix(c(1:12), nrow=3)
#' expit(exmat)
#'
#'

expit <- function(x){
  return(1/(1+exp(-x)))
}
