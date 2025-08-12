#' Deprecated. Please use plot_pairvar().
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param color Color to use in plotting.
#' @param xvar Character string for the name of the trait to be used. Must match the trait name in
#'    the input data used to build the model.
#' @param xlabel Optional character string to replace xvar when plotting.
#' @param orig.scale Logical value indicating whether to back-transform trait data to the original
#'    scale (TRUE) or leave them with mean zero and unit variance (FALSE).
#' @param ci_width A real number (0,1) describing the desired widths of credible band. Defaults to 0.95.
#' @param ymin Real number indicating the location of the bottom of the plot's y axis.
#' @param ymax Real number indicating the location of the top of the plot's y axis.
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
#' scatter_pairvar(ex_compnet_phylo, xvar="phylodist", ymax=0.25)

scatter_pairvar <- function(mod,
                            xvar,
                            xlabel,
                            color="red",
                            orig.scale=TRUE,
                            ymin=0,
                            ymax=1,
                            ci_width=0.95,
                            grid_size=100,
                            thin=TRUE,
                            thin_to=100){
  stop("Deprecated. Please use plot_pairvar().")
}
