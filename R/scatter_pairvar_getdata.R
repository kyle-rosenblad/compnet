#' Deprecated. Please use plot_pairvar_getdata().
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param xvar Character string for the name of the trait to be used. Must match the trait name in
#'    the input data used to build the model.
#' @param orig.scale Logical value indicating whether to back-transform trait data to the original
#'    scale (TRUE) or leave them with mean zero and unit variance (FALSE).
#' @param ci_width A real number (0,1) describing the desired widths of credible band. Defaults to 0.95.
#' @param grid_size A positive integer defining the number of discrete steps to use in approximating
#'    the shape of mean prediction curves and credible bands. Defaults to 100.
#' @param thin Logical value determining whether to use a random subsample of the full posterior sample.
#' @param thin_to Integer value determining how many random samples to draw from the full posterior sample if thin=TRUE.
#' @return A ggplot2 graphic.
#' @examples
#' data(ex_presabs)
#' data(ex_phylo)
#'
#' ex_compnet_phylo <- buildcompnet(presabs=ex_presabs, pairvars=ex_phylo, warmup=10, iter=20)
#'
#' scatter_pairvar_getdata(ex_compnet_phylo, xvar="phylodist")

scatter_pairvar_getdata <- function(){
  stop("Deprecated. Please use plot_pairvar_getdata().")
}
