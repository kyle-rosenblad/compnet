#' Compute summary statistics quantifying row/column-level and third order dependencies in a
#'    symmetric matrix with NA diagonals
#'
#' @importFrom stats sd
#' @export
#' @param Y A symmetric matrix with NA diagonals
#' @return A named vector containing: 1- the standard deviation of row means, and 2- The triadic
#'    dependency metric used by Hoff, Fosdick, & Volfovsky's "amen" package.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- compnet(presabs=ex_presabs, spvars_dist_int=ex_traits, warmup=100, iter=200)
#'
#' ex_compnet_summ <- summarize_compnet(ex_compnet)
#' ex_compnet_summ
#'

compute_gofstats<-function(Y)
{
  sd.rowmean <- stats::sd(rowMeans(Y,na.rm=TRUE) ,na.rm=TRUE)

  E<-Y-mean(Y,na.rm=TRUE)
  D<-1*(!is.na(E)) ; E[is.na(E)]<-0
  triad.dep<- sum(diag(E%*%E%*%E))/( sum(diag(D%*%D%*%D))*stats::sd(c(Y),na.rm=TRUE)^3)

  gof<-c(sd.rowmean, triad.dep )
  gof[is.na(gof)]<-0
  names(gof)<-c("sd.rowmean","triad.dep")
  return(gof)
}
