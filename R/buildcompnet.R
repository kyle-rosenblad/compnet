#' Network models of interspecific competition with presence-absence data
#'
#' @importFrom utils combn
#' @importFrom stats sd
#' @export
#' @param presabs Must be specified by user. Binary (0 or 1) presence-absence matrix with sites as
#'    rows and species as columns. Column names should be unique species names.
#' @param spvars_no_int A matrix or data frame, in which rows are species and columns are traits to
#'    be included in the model as additive species effects only, with no interaction term. Row names
#'    should be unique species names, and column names should be unique trait names.
#' @param spvars_dist_int A matrix or data frame, in which rows are species and columns are traits
#'    to be included in the model with a distance interaction term (i.e., the absolute value of the
#'    difference in trait values for each species pair). Row names should be unique species names,
#'    and column names should be unique trait names.
#' @param spvars_multi_int A matrix or data frame, in which rows are species and columns are traits
#'    to be included in the model with a multiplicative interaction term (i.e., the product of the
#'    trait values for each species pair). Row names should be unique species names, and column
#'    names should be unique trait names.
#' @param spvars_cat_no_int A matrix or data frame, in which rows are species and columns are categorical traits to
#'    be included in the model as additive species effects only, with no interaction term. Row names
#'    should be unique species names, and column names should be unique trait names.
#' @param spvars_cat_int A matrix or data frame, in which rows are species and columns are categorical traits to
#'    be included in the model with a binary "same or different" interaction term. Row names
#'    should be unique species names, and column names should be unique trait names.
#' @param pairvars A matrix or data frame, in which rows are species pairs and columns are
#'    pair-level traits that do not have a species-level analog, e.g., phylogenetic distance.
#'    Column names should be unique trait names. There
#'    should also be two columns named "spAid" and "spBid" containing the unique names of the
#'    species in each pair.
#' @param rank Number of dimensions for the multiplicative latent factor term. Rank=0 (the default)
#'    yields a model with no multiplicative term.
#' @param prior_intercept_scale Scale parameter for mean-zero Gaussian prior on the intercept term
#'    for the linear predictor.
#' @param prior_betas_scale Scale parameter for mean-zero Gaussian priors on the coefficients of
#'    fixed effect terms in the linear predictor.
#' @param prior_sigma_addeff_rate Rate parameter for exponential prior on the scale of the
#'    species-level Gaussian random effects (i.e., "row and column effects").
#' @param prior_lambda_scale Scale parameter for mean-zero Gaussian prior on diagonal values of
#'    Lambda, the matrix that determines how different species' values of the latent factors
#'    interact in the linear predictor.
#' @param warmup Number of warmup iterations for Stan.
#' @param iter Number of posterior sampling iterations for Stan.
#' @param adapt_delta A parameter that tunes Stan's posterior sampling algorithm. Increasing closer to 1 can help avoid divergent transitions.
#' @return Object of class "compnet", which is a list containing the stanfit model object,
#'    a named list of posterior samples for all model parameters, a data frame containing all
#'    input variables for the model, a matrix of dyadic X variables, a matrix of X variables
#'    pertaining to species A in each pair, a matrix of X variables pertaining to species B in
#'    each pair, and--when relevant--a matrix
#'    of means and standard deviations for the input trait data before centering and scaling.
#' @details This function uses Stan, as accessed through the rstan package, to build a network
#'    regression model of interspecific competitive niche differentiation in a Bayesian framework.
#'    This function is designed to test the hypothesis that species are more likely to co-occur with
#'    other species that are functionally dissimilar. Functional dissimilarity can be represented
#'    directly by traits or by a proxy like phylogenetic distance. The core input data are a
#'    species-by-site presence-absence matrix and one or more species-level traits
#'    (e.g., plant leaf size) or pair-level traits (e.g., phylogenetic distance). Units of analysis
#'    are species pairs. The response variable follows a binomial
#'    distribution. The number of trials is the number of sites occupied by at least one
#'    species in a pair, and the number of successes is the number of sites occupied by both species.
#'    If species-level traits are used, each trait can be non-interacting
#'    (i.e., there is no interaction term between species A's trait value and species B's),
#'    interacting via a typical multiplicative term, or interacting via an absolute value difference
#'    (i.e., "distance") term. The interaction terms are key to the core hypothesis. If competitive
#'    niche differentiation is occurring, then the probability of co-occurrence is expected to
#'    increase with trait or phylogenetic distance between species A and B, or with the product of
#'    their trait values, if a multiplicative interaction is specified instead of a distance interaction.
#'    Random effects are used to account for additive species-level dependencies and, optionally,
#'    higher-order dependencies involving multiple species (e.g., "the enemy of my enemy is my friend").
#'    Higher-order dependencies are modeled using a number of latent variable specified by "rank".
#'    For more details on the random effects, see Hoff, P. (2021) Additive and multiplicative effects
#'    network models. Stat. Sci. 36, 34–50.
#' @examples
#'
#' data(ex_presabs)
#' data(ex_traits)
#'
#' # Quick demo run. Will prompt warnings.
#' # Run with default warmup and iter for good posterior sampling.
#' ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
#'

buildcompnet <- function(presabs,
                    spvars_no_int=NULL,
                    spvars_dist_int=NULL,
                    spvars_multi_int=NULL,
                    spvars_cat_no_int=NULL,
                    spvars_cat_int=NULL,
                    pairvars=NULL,
                    rank=0,
                    prior_intercept_scale=5,
                    prior_betas_scale=5,
                    prior_sigma_addeff_rate=1,
                    prior_lambda_scale=5,
                    warmup=1000,
                    iter=2000,
                    adapt_delta=0.99
){

  ### errors for incorrectly formatted input data
  if(!missing("spvars_no_int")){
    if(!identical(sort(rownames(spvars_no_int)), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }

  if(!missing("spvars_multi_int")){
    if(!identical(sort(rownames(spvars_multi_int)), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }
  if(!missing("spvars_dist_int")){
    if(!identical(sort(rownames(spvars_dist_int)), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }

  if(!missing("spvars_cat_no_int")){
    if(!identical(sort(rownames(spvars_cat_no_int)), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }

  if(!missing("spvars_cat_int")){
    if(!identical(sort(rownames(spvars_cat_int)), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }

  if(!missing("pairvars")){
    tmpsplist <- unique(c(as.data.frame(pairvars)$spAid, as.data.frame(pairvars)$spBid))
    if(!identical(sort(tmpsplist), sort(colnames(presabs)))){
      stop("Please check the formatting of your presence-absence and trait data. see ?buildcompnet()")
    }
  }

  ### data prep
  dyads <- t(utils::combn(colnames(presabs),2))
  d <- as.data.frame(dyads)
  names(d) <- c("spAid", "spBid")

  for(j in 1:nrow(d)){
    tmp <- presabs[,dyads[j,c(1:2)]]
    d[j, "both"] <- length(which(tmp[,1]==1 & tmp[,2]==1))
    d[j, "either"] <- length(which(tmp[,1]==1 | tmp[,2]==1))
  }

  XA <- matrix(NA, ncol=1, nrow=nrow(d))
  XB <- matrix(NA, ncol=1, nrow=nrow(d))
  Xdy <- matrix(NA, ncol=1, nrow=nrow(d))

  if(!missing("spvars_no_int")){
    for(i in 1:ncol(spvars_no_int)){
      tmptrait <- scale(spvars_no_int[,i])
      names(tmptrait) <- rownames(spvars_no_int)
      vecA <- c()
      vecB <- c()
      for(j in 1:nrow(d)){
        vecA[j] <- tmptrait[d[j, "spAid"]]
        vecB[j] <- tmptrait[d[j, "spBid"]]
      }
      XA <- as.matrix(cbind(XA, vecA))
      XB <- as.matrix(cbind(XB, vecB))
      colnames(XA)[length(colnames(XA))] <- paste(names(spvars_no_int[i]), "A", sep="_")
      colnames(XB)[length(colnames(XB))] <- paste(names(spvars_no_int[i]), "B", sep="_")
    }
  }

  if(!missing("spvars_dist_int")){
    for(i in 1:ncol(spvars_dist_int)){
      tmptrait <- scale(spvars_dist_int[,i])
      names(tmptrait) <- rownames(spvars_dist_int)
      vecA <- c()
      vecB <- c()
      vecdy <- c()
      for(j in 1:nrow(d)){
        vecA[j] <- tmptrait[d[j, "spAid"]]
        vecB[j] <- tmptrait[d[j, "spBid"]]
        vecdy[j] <- abs(tmptrait[d[j, "spBid"]]-tmptrait[d[j, "spAid"]])
      }
      XA <- as.matrix(cbind(XA, vecA))
      XB <- as.matrix(cbind(XB, vecB))
      Xdy <- as.matrix(cbind(Xdy, vecdy))
      colnames(XA)[length(colnames(XA))] <- paste(names(spvars_dist_int[i]), "A", sep="_")
      colnames(XB)[length(colnames(XB))] <- paste(names(spvars_dist_int[i]), "B", sep="_")
      colnames(Xdy)[length(colnames(Xdy))] <- paste(names(spvars_dist_int[i]), "dist", sep="_")
    }
  }

  if(!missing("spvars_multi_int")){
    for(i in 1:ncol(spvars_multi_int)){
      tmptrait <- scale(spvars_multi_int[,i])
      names(tmptrait) <- rownames(spvars_multi_int)
      vecA <- c()
      vecB <- c()
      vecdy <- c()
      for(j in 1:nrow(d)){
        vecA[j] <- tmptrait[d[j, "spAid"]]
        vecB[j] <- tmptrait[d[j, "spBid"]]
        vecdy[j] <- tmptrait[d[j, "spBid"]]*tmptrait[d[j, "spAid"]]
      }
      XA <- as.matrix(cbind(XA, vecA))
      XB <- as.matrix(cbind(XB, vecB))
      Xdy <- as.matrix(cbind(Xdy, vecdy))
      colnames(XA)[length(colnames(XA))] <- paste(names(spvars_multi_int[i]), "A", sep="_")
      colnames(XB)[length(colnames(XB))] <- paste(names(spvars_multi_int[i]), "B", sep="_")
      colnames(Xdy)[length(colnames(Xdy))] <- paste(names(spvars_multi_int[i]), "multi", sep="_")
    }
  }

  if(!missing("spvars_cat_no_int")){
    for(i in 1:ncol(spvars_cat_no_int)){
      tmptrait <- factor(spvars_cat_no_int[,i])
      names(tmptrait) <- rownames(spvars_cat_no_int)
      cats <- levels(tmptrait)
      dummiesA <- matrix(NA, nrow=nrow(d), ncol=length(cats))
      dummiesB <- matrix(NA, nrow=nrow(d), ncol=length(cats))
      for(k in 1:length(cats)){
        for(j in 1:nrow(d)){
          dummiesA[j, k] <- as.numeric(tmptrait[d[j, "spAid"]] == cats[k])
          dummiesB[j, k] <- as.numeric(tmptrait[d[j, "spBid"]] == cats[k])
        }
      }
      colnames(dummiesA) <- paste(names(spvars_cat_no_int[i]), cats, "dummy_A", sep="_")
      colnames(dummiesB) <- paste(names(spvars_cat_no_int[i]), cats, "dummy_B", sep="_")
      XA <- as.matrix(cbind(XA, dummiesA[,-1])) # drop first category as reference
      XB <- as.matrix(cbind(XB, dummiesB[,-1])) # drop first category as reference
    }
  }

  if(!missing("spvars_cat_int")){
    for(i in 1:ncol(spvars_cat_int)){
      tmptrait <- factor(spvars_cat_int[,i])
      names(tmptrait) <- rownames(spvars_cat_int)
      cats <- levels(tmptrait)
      dummiesA <- matrix(NA, nrow=nrow(d), ncol=length(cats))
      dummiesB <- matrix(NA, nrow=nrow(d), ncol=length(cats))
      for(k in 1:length(cats)){
        for(j in 1:nrow(d)){
          dummiesA[j, k] <- as.numeric(tmptrait[d[j, "spAid"]] == cats[k])
          dummiesB[j, k] <- as.numeric(tmptrait[d[j, "spBid"]] == cats[k])
        }
      }
      vecdy <- c()
      for(j in 1:nrow(d)){
        vecdy[j] <- as.numeric(tmptrait[d[j, "spBid"]] == tmptrait[d[j, "spAid"]])
      }
      colnames(dummiesA) <- paste(names(spvars_cat_int[i]), cats, "dummy_A", sep="_")
      colnames(dummiesB) <- paste(names(spvars_cat_int[i]), cats, "dummy_B", sep="_")
      XA <- as.matrix(cbind(XA, dummiesA[,-1])) # drop first category as reference
      XB <- as.matrix(cbind(XB, dummiesB[,-1])) # drop first category as reference
      Xdy <- as.matrix(cbind(Xdy, vecdy))
      colnames(Xdy)[length(colnames(Xdy))] <- paste(names(spvars_cat_int[i]), "int", sep="_")
    }
  }

  if(!missing("pairvars")){
    pairvars.orig <- as.data.frame(pairvars)
    pairvars.orig[c("spAid", "spBid")] <- NULL
    pairvars <- as.matrix(pairvars)
    for(j in 1:nrow(pairvars)){
      pairvars[j, c("spAid", "spBid")] <- sort(pairvars[j, c("spAid", "spBid")])
    }
    pairvars <- as.data.frame(pairvars)
    pairvars <- pairvars[order(pairvars$spBid),]
    pairvars <- pairvars[order(pairvars$spAid),]
    pairvars$spAid <- NULL
    pairvars$spBid <- NULL
    for(k in 1:ncol(pairvars)){
      pairvars[,k] <- scale(as.numeric(pairvars[,k]))
    }
    Xdy <- as.matrix(cbind(Xdy, pairvars))
  }

  ### deal with column names in X matrices for special case where we only have one X variable
  if(ncol(Xdy)==2){
    Xdyname <- colnames(Xdy)[2]
  }
  if(ncol(XA)==2){
    XAname <- colnames(XA)[2]
  }
  if(ncol(XB)==2){
    XBname <- colnames(XB)[2]
  }

  Xdy <- as.matrix(Xdy[,-1])
  XA <- as.matrix(XA[,-1])
  XB <- as.matrix(XB[,-1])

  if(ncol(Xdy)==1){
    colnames(Xdy) <- Xdyname
  }
  if(ncol(XA)==1){
    colnames(XA) <- XAname
  }
  if(ncol(XB)==1){
    colnames(XB) <- XBname
  }
  # save the summary statistics needed to return continuous x variables to original scales
  if(!missing(spvars_no_int)){
    spvars_no_int_summs <- as.data.frame(cbind(apply(spvars_no_int, 2, mean), apply(spvars_no_int, 2, stats::sd)))
    names(spvars_no_int_summs) <- c("mean", "sd")
  }

  if(!missing(spvars_dist_int)){
    spvars_dist_summs <- as.data.frame(cbind(apply(spvars_dist_int, 2, mean), apply(spvars_dist_int, 2, stats::sd)))
    names(spvars_dist_summs) <- c("mean", "sd")
  }

  if(!missing(spvars_multi_int)){
    spvars_multi_summs <- as.data.frame(cbind(apply(spvars_multi_int, 2, mean), apply(spvars_multi_int, 2, stats::sd)))
    names(spvars_multi_summs) <- c("mean", "sd")
  }

  if(!missing(pairvars)){
    pairvars_summs <- as.data.frame(cbind(apply(pairvars.orig, 2, mean), apply(pairvars.orig, 2, stats::sd)))
    names(pairvars_summs) <- c("mean", "sd")
  }

  # convert species IDs to integers for Stan
  d$spAid_orig <- d$spAid
  d$spBid_orig <- d$spBid
  d$spAid <- match(d$spAid, colnames(presabs))
  d$spBid <- match(d$spBid, colnames(presabs))

  # save data frame of model inputs
  input_df <- as.data.frame(cbind(d[c("spAid", "spBid", "spAid_orig", "spBid_orig", "both", "either")], Xdy, XA, XB))



  ### build model
  if(rank==0){
    datalist <- list(
      n_nodes=ncol(presabs),
      N=nrow(d),
      Xdy_cols=ncol(Xdy),
      Xsp_cols=ncol(XA),
      spAid=d$spAid,
      spBid=d$spBid,
      Xdy=Xdy,
      XA=XA,
      XB=XB,
      both=d$both,
      either=d$either,
      prior_intercept_scale=prior_intercept_scale,
      prior_betas_scale=prior_betas_scale,
      prior_sigma_addeff_rate=prior_sigma_addeff_rate)

    stanmod <- rstan::sampling(stanmodels$srm_binomial,
                               data=datalist,
                               cores=1,
                               chains=1,
                               warmup=warmup,
                               iter=iter,
                               verbose=F,
                               control=list(adapt_delta=adapt_delta))
  }

  if(rank>0){
    datalist <- list(
      n_nodes=ncol(presabs),
      N=nrow(d),
      Xdy_cols=ncol(Xdy),
      Xsp_cols=ncol(XA),
      spAid=d$spAid,
      spBid=d$spBid,
      Xdy=Xdy,
      XA=XA,
      XB=XB,
      K=rank,
      both=d$both,
      either=d$either,
      prior_intercept_scale=prior_intercept_scale,
      prior_betas_scale=prior_betas_scale,
      prior_sigma_addeff_rate=prior_sigma_addeff_rate,
      prior_lambda_scale=prior_lambda_scale)

    stanmod <- rstan::sampling(stanmodels$ame_binomial,
                               data=datalist,
                               pars=c("intercept", # Marginalize over the individual values of each latent factor ("U"), since these are known to be nonidentifiable due to reflectional and rotational invariance, as well as potential label-switching among multiple latent factor dimensions when K>1. The Stan code has tools designed to deal with this, but only for computational efficiency reasons; it's not a problem for inference.
                                      "beta_dy",
                                      "beta_sp",
                                      "sigma_species_randint",
                                      "species_randint",
                                      "xb",
                                      "pboth",
                                      "latfacterm",
                                      "lp__"),
                               cores=1,
                               chains=1,
                               warmup=warmup,
                               iter=iter,
                               verbose=F,
                               control=list(adapt_delta=adapt_delta))
  }



  ### extract posterior samples
  stanmod_samp <- rstan::extract(stanmod)


  ### assemble list of outputs
  outlist <- list(stanmod, stanmod_samp, input_df, Xdy, XA, XB)
  names(outlist) <- c("stanmod", "stanmod_samp", "d", "Xdy", "XA", "XB")

  if(!missing(spvars_no_int)){
    outlist[[length(outlist)+1]] <- spvars_no_int_summs
    names(outlist)[length(outlist)] <- "spvars_no_int_summs"
  }

  if(!missing(spvars_dist_int)){
    outlist[[length(outlist)+1]] <- spvars_dist_summs
    names(outlist)[length(outlist)] <- "spvars_dist_summs"
  }

  if(!missing(spvars_multi_int)){
    outlist[[length(outlist)+1]] <- spvars_multi_summs
    names(outlist)[length(outlist)] <- "spvars_multi_summs"
  }

  if(!missing(pairvars)){
    outlist[[length(outlist)+1]] <- pairvars_summs
    names(outlist)[length(outlist)] <- "pairvars_summs"
  }

  print(paste0("compnet uses Stan under the hood. When Stan finishes fitting your model, ",
        "it may issue warnings. To deal with any warnings Stan might issue, ",
        "Please see the links provided in Stan's output, as well as the compnet website:",
        "https://kyle-rosenblad.github.io/compnet/"))

  class(outlist) <- "compnet"
  return(outlist)
}
