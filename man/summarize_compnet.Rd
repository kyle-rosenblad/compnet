% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_compnet.R
\name{summarize_compnet}
\alias{summarize_compnet}
\title{Summarizing compnet model output}
\usage{
summarize_compnet(mod, ci_width = 0.95)
}
\arguments{
\item{mod}{Object of class "compnet", which is created by the buildcompnet() function.}

\item{ci_width}{A real number (0,1) of the desired interval width. Defaults to 0.95.}
}
\value{
A data frame summarizing means and credible intervals for standardized effect sizes of fixed effects.
}
\description{
Summarizing compnet model output
}
\examples{

data(ex_presabs)
data(ex_traits)

# Quick demo run. Will prompt warnings.
# Run with default warmup and iter for good posterior sampling.
ex_compnet <- buildcompnet(presabs=ex_presabs, spvars_dist_int=ex_traits[c("ndtrait")], warmup=10, iter=20)

ex_compnet_summ <- summarize_compnet(ex_compnet)
ex_compnet_summ

}
