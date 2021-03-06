% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateAUCwithPairwiseCML.R
\encoding{UTF-8}
\name{estimateAUCwithPairwiseCML}
\alias{estimateAUCwithPairwiseCML}
\title{estimate AUCwith pairwise censored maximum likelihood}
\usage{
estimateAUCwithPairwiseCML(
  inputData,
  LOQ,
  timePoints,
  isMultiplicative = FALSE,
  onlyFitCML = FALSE,
  optimizationMethod = NULL,
  CMLcontrol = NULL,
  na.rm = TRUE
)
}
\arguments{
\item{inputData}{numeric matrix or data frame of the size 
n by J (n the sample size and J the number of time points)
 the input dataset}

\item{LOQ}{scalar, limit of quantification value}

\item{timePoints}{vector of time points}

\item{isMultiplicative}{logical variable indicating whether
an additive error model (FALSE) or a multiplicative error
model (TRUE) should be used}

\item{onlyFitCML}{logical variable with FALSE as default, if TRUE only the
censored maximum likelihood estimates will be calculated.}

\item{optimizationMethod}{single string specifying the method to be used for optimizing the log-likelihood, 
the default is NULL that allows the function to decide the about the best method. Otherwise, one can select among choices
available via R package maxLik: "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), 
"BFGSR" (for the BFGS algorithm implemented in R), 
"BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), 
"CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). 
Lower-case letters (such as "nr" for Newton-Raphson) are allowed.}

\item{CMLcontrol}{list of arguments to control 
convergence of maximization algorithm. It is the same argument
as control in the function maxLik in the R package maxLik}

\item{na.rm}{logical variable indicating whether the lines with missing values
should be ignored (TRUE, default) or not (FALSE). Note that, it will be applied 
for the sub-datasets regarding each pair.}
}
\value{
a list with three components: output of maxLik function,
estimated parameters (mean vector and the covariance matrix) 
using censored maximum likelihood, and estimated
AUC and its standard error.
}
\description{
function to estimate mean and and covariance matrix of censored data using a 
full censored maximum likelihood approach via fitting all possible pairs, then use these 
estimates for estimating AUC and its standard error
}
\examples{
# generate data from Beal model with only fixed effects
set.seed(111)
genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
		1, 1, seq(0.5,3,1.5))
estimateAUCwithPairwiseCML(genDataFixedEffects, 0.1, seq(0.5,3,1.5))
}
\seealso{
\href{https://www.rdocumentation.org/packages/maxLik/versions/1.3-4/topics/maxLik}{maxLik}
}
\author{
Vahid Nassiri, Helen Yvette Barnett
}
