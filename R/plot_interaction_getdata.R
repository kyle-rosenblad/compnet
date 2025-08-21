#' Get data for a plot of the interactive effect of two species' values of the same trait.
#'
#' @importFrom stats quantile
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
#' @return A data frame containing data needed to plot mean curves and credible bands as done in
#'    the plot_interaction function.
#' @examples
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- buildcompnet(presabs=ex_presabs,
#' spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20, family='binomial')
#'
#' plotdata <- plot_interaction_getdata(ex_compnet, xvar="ndtrait")

plot_interaction_getdata <- function(mod,
                                     xvar,
                                     orig.scale=T,
                                     intlevels=c(0.05,0.5,0.95),
                                     ci_width=0.95,
                                     grid_size=100,
                                     thin=TRUE,
                                     thin_to=100){
  ### error for categorical xvar or pairvar or non-interacting continuous trait
  if(xvar%in%rownames(mod$spvars_multi_summs)==FALSE &
     xvar%in%rownames(mod$spvars_dist_summs)==FALSE){
    stop("This function only supports continuous traits with multiplicative or distance interactions. see ?plot_interaction_getdata()")
  }

  intlevels <- sort(intlevels)

  lowquant <- (1-ci_width)/2
  highquant <- 1-lowquant
  samp_temp <- mod$stanmod_samp

  # for the trait of interest get posterior samples of:
  # intercept, species beta, and interaction beta (either distance or product).
  alpha <- samp_temp$intercept
  colpos <- which(colnames(mod$XA)==paste(xvar, "A", sep="_"))
  beta_sp <- samp_temp$beta_sp[,colpos]

  # get other betas:
  beta_sp_other <- as.matrix(samp_temp$beta_sp[,-colpos])

  # get x data for trait of interest
  XAmin <- min(mod$XA[,colpos])
  XAmax <- max(mod$XA[,colpos])
  XBmin <- min(mod$XB[,colpos])
  XBmax <- max(mod$XB[,colpos])

  # get x values for interaction plot levels
  intlevels <- XBmin+intlevels*(XBmax-XBmin)

  # get mean values of other x variables:
  XA_other_mean <- apply(as.matrix(mod$XA[,-colpos]), 2, mean)
  XB_other_mean <- apply(as.matrix(mod$XB[,-colpos]), 2, mean)

  # get beta and x values for dyadic variables
  if(xvar%in%rownames(mod$spvars_dist_summs)){
    colpos <- which(colnames(mod$Xdy)==paste(xvar, "dist", sep="_"))
    beta_dy <- samp_temp$beta_dy[,colpos]
    beta_dy_other <- as.matrix(samp_temp$beta_dy[,-colpos])
    Xdy_other_mean <- as.matrix(apply(as.matrix(mod$Xdy[,-colpos]), 2, mean))
  }

  if(xvar%in%rownames(mod$spvars_multi_summs)){
    colpos <- which(colnames(mod$Xdy)==paste(xvar, "multi", sep="_"))
    beta_dy <- samp_temp$beta_dy[,colpos]
    beta_dy_other <- as.matrix(samp_temp$beta_dy[,-colpos])
    Xdy_other_mean <- as.matrix(apply(as.matrix(mod$Xdy[,-colpos]), 2, mean))
  }

  xbeta_other_A <- beta_sp_other%*%XA_other_mean
  xbeta_other_B <- beta_sp_other%*%XB_other_mean
  xbeta_other_dy <- beta_dy_other%*%Xdy_other_mean

  xbeta_other <- xbeta_other_A + xbeta_other_B + xbeta_other_dy

  samp_for_plot <- data.frame(alpha=alpha, beta_sp=beta_sp, beta_dy=beta_dy, xbeta_other=xbeta_other)
  grid_for_plot <- data.frame(
    x=seq(
      from=XAmin,
      to=XAmax,
      length.out=grid_size))

  if(thin==TRUE){
    samp_for_plot <- samp_for_plot[sample(1:nrow(samp_for_plot), replace=FALSE, size=min(nrow(samp_for_plot),thin_to)),]
  }



  if(mod$family=='binomial'){
    for(k in 1:length(intlevels)){
      gridtmp <- grid_for_plot
      for(i in 1:nrow(grid_for_plot)){
        x_tmp <- grid_for_plot[i, "x"]
        yvec <- c()
        for(j in 1:nrow(samp_for_plot)){
          intercept_tmp <- samp_for_plot[j, "intercept"]
          beta_sp_tmp <- samp_for_plot[j, "beta_sp"]
          beta_dy_tmp <- samp_for_plot[j, "beta_dy"]
          xbeta_other_tmp <- samp_for_plot[j, "xbeta_other"]

          if(xvar%in%rownames(mod$spvars_dist_summs)){
            ytmp <- expit(intercept_tmp +
                            xbeta_other_tmp +
                            beta_sp_tmp*intlevels[k] +
                            beta_sp_tmp*grid_for_plot[i,"x"] +
                            beta_dy_tmp*abs(grid_for_plot[i,"x"] - intlevels[k]))
          }
          if(xvar%in%rownames(mod$spvars_multi_summs)){
            ytmp <- expit(intercept_tmp +
                            xbeta_other_tmp +
                            beta_sp_tmp*intlevels[k] +
                            beta_sp_tmp*grid_for_plot[i,"x"] +
                            beta_dy_tmp*grid_for_plot[i,"x"]*intlevels[k])
          }

          yvec <- c(yvec, ytmp)
        }
        gridtmp[i, "qlow"] <- stats::quantile(yvec, lowquant)
        gridtmp[i, "means"] <- mean(yvec)
        gridtmp[i, "qhigh"] <- stats::quantile(yvec, highquant)
      }
      gridtmp$intlevel <- intlevels[k]

      if(k==1){
        gridfinal <- gridtmp
      }
      if(k>1){
        gridfinal <- rbind(gridfinal, gridtmp)
      }
    }
  }


  if(mod$family=='fnchypg'){
    for(k in 1:length(intlevels)){
      gridtmp <- grid_for_plot
      for(i in 1:nrow(grid_for_plot)){
        x_tmp <- grid_for_plot[i, "x"]
        yvec <- c()
        for(j in 1:nrow(samp_for_plot)){
          intercept_tmp <- samp_for_plot[j, "intercept"]
          beta_sp_tmp <- samp_for_plot[j, "beta_sp"]
          beta_dy_tmp <- samp_for_plot[j, "beta_dy"]
          xbeta_other_tmp <- samp_for_plot[j, "xbeta_other"]

          if(xvar%in%rownames(mod$spvars_dist_summs)){
            ytmp <- (intercept_tmp +
                       xbeta_other_tmp +
                       beta_sp_tmp*intlevels[k] +
                       beta_sp_tmp*grid_for_plot[i,"x"] +
                       beta_dy_tmp*abs(grid_for_plot[i,"x"] - intlevels[k]))
          }
          if(xvar%in%rownames(mod$spvars_multi_summs)){
            ytmp <- (intercept_tmp +
                       xbeta_other_tmp +
                       beta_sp_tmp*intlevels[k] +
                       beta_sp_tmp*grid_for_plot[i,"x"] +
                       beta_dy_tmp*grid_for_plot[i,"x"]*intlevels[k])
          }

          yvec <- c(yvec, ytmp)
        }
        gridtmp[i, "qlow"] <- stats::quantile(yvec, lowquant)
        gridtmp[i, "means"] <- mean(yvec)
        gridtmp[i, "qhigh"] <- stats::quantile(yvec, highquant)
      }
      gridtmp$intlevel <- intlevels[k]

      if(k==1){
        gridfinal <- gridtmp
      }
      if(k>1){
        gridfinal <- rbind(gridfinal, gridtmp)
      }
    }
  }


  d <- mod$d
  if(orig.scale==T){
    if(xvar%in%rownames(mod$spvars_dist_summs)){
      gridfinal$x <- gridfinal$x*mod$spvars_dist_summs[xvar,"sd"] + mod$spvars_dist_summs[xvar,"mean"]
      gridfinal$intlevel <- gridfinal$intlevel*mod$spvars_dist_summs[xvar,"sd"] + mod$spvars_dist_summs[xvar,"mean"]
      d[[paste(xvar, "A", sep="_")]] <- d[[paste(xvar, "A", sep="_")]]*mod$spvars_dist_summs[xvar,"sd"] + mod$spvars_dist_summs[xvar,"mean"]
    }
    if(xvar%in%rownames(mod$spvars_multi_summs)){
      gridfinal$x <- gridfinal$x*mod$spvars_multi_summs[xvar,"sd"] + mod$spvars_multi_summs[xvar,"mean"]
      gridfinal$intlevel <- gridfinal$intlevel*mod$spvars_multi_summs[xvar,"sd"] + mod$spvars_multi_summs[xvar,"mean"]
      d[[paste(xvar, "A", sep="_")]] <- d[[paste(xvar, "A", sep="_")]]*mod$spvars_multi_summs[xvar,"sd"] + mod$spvars_multi_summs[xvar,"mean"]
    }
  }

  return(gridfinal)
}
