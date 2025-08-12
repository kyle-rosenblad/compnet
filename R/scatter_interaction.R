#' Deprecated. Please use plot_interaction().
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param xvar Character string for the name of the trait to be used. Must match the trait name in
#'    the input data used to build the model.
#' @param xlabel Optional character string to replace xvar when plotting.
#' @param orig.scale Logical value indicating whether to back-transform trait data to the original
#'    scale (TRUE) or leave them with mean zero and unit variance (FALSE).
#' @param intlevels Vector of real values on the interval \eqn{[0,1]} indicating what levels of the x
#'    variable to condition on for species B when plotting species A's mean response.
#' @param ci_width A real number (0,1) describing the desired widths of credible bands. Defaults to 0.95.
#' @param ymin Real number indicating the location of the bottom of the plot's y axis.
#' @param ymax Real number indicating the location of the top of the plot's y axis.
#' @param grid_size A positive integer defining the number of discrete steps to use in approximating
#'    the shape of mean prediction curves and credible bands. Defaults to 100.
#' @param thin Logical value determining whether to use a random subsample of the full posterior sample.
#' @param thin_to Integer value determining how many random samples to draw from the full posterior sample if thin=TRUE.
#' @return A ggplot2 graphic.
#' @examples
#' data(ex_presabs)
#' data(ex_traits)
#'
#' ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#' plotdata <- scatter_interaction(ex_compnet, xvar="ndtrait")

#library(ggplot2)
scatter_interaction <- function(){
  stop("Deprecated. Please use plot_interaction().")
}
