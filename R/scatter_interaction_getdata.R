#' Deprecated. Please use plot_interaction_getdata().
#'
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @param xvar Character string for the name of the trait to be used. Must match the trait name in
#'    the input data used to build the model.
#' @param orig.scale Logical value indicating whether to back-transform trait data to the original
#'    scale (T) or leave them with mean zero and unit variance (F).
#' @param intlevels Vector of real values on the interval \eqn{[0,1]} indicating what levels of the x
#'    variable to condition on for species B when plotting species A's mean response.
#' @param ci_width A real number (0,1) describing the desired widths of credible bands. Defaults to 0.95.
#' @param grid_size A positive integer defining the number of discrete steps to use in approximating
#'    the shape of mean prediction curves and credible bands. Defaults to 100.
#' @param thin Logical value determining whether to use a random subsample of the full posterior sample.
#' @param thin_to Integer value determining how many random samples to draw from the full posterior
#'    sample if thin=TRUE.

scatter_interaction_getdata <- function(mod,
                                        xvar,
                                        orig.scale=T,
                                        intlevels=c(0.05,0.5,0.95),
                                        ci_width=0.95,
                                        grid_size=100,
                                        thin=TRUE,
                                        thin_to=100){
  stop("Deprecated. Please use plot_interaction_getdata().")
}
