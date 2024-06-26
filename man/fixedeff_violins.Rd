% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixedeff_violins.R
\name{fixedeff_violins}
\alias{fixedeff_violins}
\title{Violin plots of fixed effect coefficients in a compnet model}
\usage{
fixedeff_violins(mod)
}
\arguments{
\item{mod}{An object of class "compnet" created by the buildcompnet() function.}
}
\value{
A ggplot2 graphic showing violin plots of standarized effect sizes for fixed effects in a compnet model.
}
\description{
Violin plots of fixed effect coefficients in a compnet model
}
\examples{

data(ex_presabs)
data(ex_traits)

# Quick demo run. Will prompt warnings.
# Run with default warmup and iter for good posterior sampling.
ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)
fixedeff_violins(ex_compnet)

}
