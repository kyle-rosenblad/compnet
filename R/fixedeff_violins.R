#' Violin plots of fixed effect coefficients in a compnet model
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_violin xlab ylab theme_bw
#' @export
#' @param mod An object of class "compnet" created by the buildcompnet() function.
#' @return A ggplot2 graphic showing violin plots of standarized effect sizes for fixed effects in a compnet model.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- buildcompnet(presabs=ex_presabs,
#' spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#' fixedeff_violins(ex_compnet)
#'

fixedeff_violins <- function(mod){
  samp <- as.data.frame(mod$stanmod_samp$intercept)
  names(samp) <- "val"
  samp$param <- "Intercept"

  if("beta_sp"%in%names(mod$stanmod_samp)){
    samp_sp_tmp <- mod$stanmod_samp$beta_sp
    for(i in 1:ncol(samp_sp_tmp)){
      if(i==1){
        samp_sp <- data.frame(param=gsub(x=colnames(mod$XA)[i], pattern="_A", replacement=""),
                              val=samp_sp_tmp[,i])
      }
      if(i>1){
        samp_sp <- rbind(samp_sp,
                         data.frame(param=gsub(x=colnames(mod$XA)[i], pattern="_A", replacement=""),
                              val=samp_sp_tmp[,i]))

      }
    }
  samp <- rbind(samp, samp_sp)
  }


  if("beta_dy"%in%names(mod$stanmod_samp)){
    samp_dy_tmp <- mod$stanmod_samp$beta_dy
    for(i in 1:ncol(samp_dy_tmp)){
      if(i==1){
        samp_dy <- data.frame(param=colnames(mod$Xdy)[i],
                              val=samp_dy_tmp[,i])
      }
      if(i>1){
        samp_dy <- rbind(samp_dy,
                         data.frame(param=colnames(mod$Xdy)[i],
                                    val=samp_dy_tmp[,i]))

      }
    }
    samp <- rbind(samp, samp_dy)
  }

  samp$param <- gsub(samp$param, pattern="_dist", replacement=" Distance")
  samp$param <- gsub(samp$param, pattern="_multi", replacement=" Product")
  samp$param <- factor(samp$param)
  samp$param <- factor(samp$param, levels=c("Intercept", levels(samp$param)[-which(levels(samp$param)=="Intercept")]))

  ggplot2::ggplot(data=samp, ggplot2::aes(x=.data$param, y=.data$val))+
    ggplot2::geom_hline(yintercept=0)+
    ggplot2::geom_violin()+
    ggplot2::xlab("Model Coefficient")+
    ggplot2::ylab("Standardized Effect Size")+
    ggplot2::theme_bw()
}
