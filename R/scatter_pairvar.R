#' Make an automated plot of the effect of a pair-level variable like phylogenetic distance
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point xlab ylab ylim theme_bw theme
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
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
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
  lowquant <- (1-ci_width)/2
  highquant <- 1-lowquant
  samp_temp <- mod$stanmod_samp

  # get posterior samples of intercept and the beta of interest
  alpha <- samp_temp$intercept
  colpos <- which(colnames(mod$Xdy)==xvar)
  beta <- samp_temp$beta_dy[,colpos]

  # get x data for trait of interest
  Xmin <- min(mod$Xdy[,colpos])
  Xmax <- max(mod$Xdy[,colpos])

  xbeta_other <- 0

  # get species-level betas and mean values of species-level variables:
  if(exists("ex_compnet_phylo$stanmod_samp$beta_sp")){
    beta_sp_other <- as.matrix(samp_temp$beta_sp)
    XA_other_mean <- apply(as.matrix(mod$XA), 2, mean)
    XB_other_mean <- apply(as.matrix(mod$XB), 2, mean)
    xbeta_other_A <- beta_sp_other%*%XA_other_mean
    xbeta_other_B <- beta_sp_other%*%XB_other_mean
    xbeta_other <- xbeta_other + xbeta_other_A + xbeta_other_B
  }

  # get beta and mean x values for other dyadic variables
  if(ncol(mod$Xdy)>1){
    beta_dy_other <- as.matrix(samp_temp$beta_dy[,-colpos])
    Xdy_other_mean <- as.matrix(apply(as.matrix(mod$Xdy[,-colpos]), 2, mean))
    xbeta_other_dy <- beta_dy_other%*%Xdy_other_mean
    xbeta_other <- xbeta_other + xbeta_other_dy
  }

  samp_for_plot <- data.frame(alpha=alpha, beta=beta, xbeta_other=xbeta_other)
  grid_for_plot <- data.frame(
    x=seq(
      from=Xmin,
      to=Xmax,
      length.out=grid_size))

  if(thin==TRUE){
    samp_for_plot <- samp_for_plot[sample(1:nrow(samp_for_plot), replace=FALSE, size=min(nrow(samp_for_plot),thin_to)),]
  }

  for(i in 1:nrow(grid_for_plot)){
    x_tmp <- grid_for_plot[i, "x"]
    yvec <- c()
    for(j in 1:nrow(samp_for_plot)){
      alpha_tmp <- samp_for_plot[j, "alpha"]
      beta_tmp <- samp_for_plot[j, "beta"]
      xbeta_other_tmp <- samp_for_plot[j, "xbeta_other"]
      ytmp <- expit(alpha_tmp +
                      xbeta_other_tmp +
                      beta_tmp*grid_for_plot[i,"x"])
      yvec <- c(yvec, ytmp)
    }
    grid_for_plot[i, "qlow"] <- stats::quantile(yvec, lowquant)
    grid_for_plot[i, "means"] <- mean(yvec)
    grid_for_plot[i, "qhigh"] <- stats::quantile(yvec, highquant)
  }
  gridfinal <- grid_for_plot

  d <- mod$d
  if(orig.scale==TRUE){
    gridfinal$x <- gridfinal$x*mod$pairvars_summs[xvar,"sd"] + mod$pairvars_summs[xvar,"mean"]
    d[[xvar]] <- d[[xvar]]*mod$pairvars_summs[xvar,"sd"] + mod$pairvars_summs[xvar,"mean"]
  }

  if(missing(xlabel)){
    xlabel <- xvar
  }

  ggplot2::ggplot()+
    ggplot2::geom_ribbon(data=gridfinal, aes(x=.data$x, ymin=.data$qlow, ymax=.data$qhigh), alpha=0.2, fill=color)+
    ggplot2::geom_line(data=gridfinal, aes(x=.data$x, y=.data$means), lwd=1, color=color)+
    ggplot2::geom_point(data=d, alpha=0.2, aes(x=.data[[xvar]], y=.data$both/.data$either))+
    ggplot2::xlab(xlabel)+
    ggplot2::ylab("P(Co-occurrence)")+
    ggplot2::ylim(c(ymin, ymax))+
    ggplot2::theme_bw()+
    ggplot2::theme(aspect.ratio=1)
}
