#' Summarizing compnet model output
#'
#' @importFrom stats quantile
#' @export
#' @param mod Object of class "compnet", which is created by the compnet function.
#' @param ci_width A real number (0,1) of the desired interval width. Defaults to 0.95.
#' @return A data frame summarizing means and credible intervals for standardized effect sizes of fixed effects.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- compnet(presabs=ex_presabs, spvars_dist_int=ex_traits, warmup=10, iter=20)
#'
#' ex_compnet_summ <- summarize_compnet(ex_compnet)
#' ex_compnet_summ
#'

summarize_compnet <- function(mod,
                              ci_width=0.95){
  lowquant <- (1-ci_width)/2
  highquant <- 1-lowquant

  samp_temp <- mod$stanmod_samp

  samp <- mod$stanmod_samp

  samp_int <- as.data.frame(samp$intercept)
  names(samp_int) <- "intercept"

  if("beta_dy"%in%names(samp)){
    samp_dy <- as.data.frame(samp$beta_dy)
    names(samp_dy) <- colnames(mod$Xdy)

  }

  if("beta_sp"%in%names(samp)){
    samp_sp <- as.data.frame(samp$beta_sp)
    names(samp_sp) <- gsub("_A", "_sp", colnames(mod$XA))
  }

  samp2 <- samp_int

  if("beta_dy"%in%names(samp)){
    samp2 <- cbind(samp2, samp_dy)
  }

  if("beta_sp"%in%names(samp)){
    samp2 <- cbind(samp2, samp_sp)
  }

  sampmeans <- as.matrix(apply(samp2, 2, mean))
  sampcis <- t(as.matrix(apply(samp2, 2, stats::quantile, probs=c(lowquant, highquant))))

  summ <- as.data.frame(cbind(sampmeans,sampcis))
  names(summ)[1] <- "Mean"
  return(summ)

}
