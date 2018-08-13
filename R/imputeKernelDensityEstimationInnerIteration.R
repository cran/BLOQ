#################################################################
### description:
###
### compuites the inner itration to produce one imputed value
### for the BLOQ observations.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' function to produce a single imputation using kernel density 
#' estimation, this will be used an in innter itreration for the
#' function imputeKernelDensityEstimation
#' @encoding UTF-8
#'  @param inputCol numeric vector of the size 
#' n (the sample size)
#' @param LOQ scalar limit of quantification value
#' @param epsilon scalar the difference between two iterations 
#' which achieving it would stop the procedure (convergence)
#' @param maxIter scalar, the maximum number of iterations with 1000 as default.
#' @return single imputed value
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @noRd
imputeKernelDensityEstimationInnerIteration <- function(
		inputCol, LOQ, epsilon, maxIter){
	## checking the input
	if (length(inputCol) < 1){
		stop("At least one non-BLOQ observation is needed to impute using kernel density estimator")
	}
	## initializing 
	imputedObs <- 0 
	## estimating the density
	estimatedDensity <- density(c(imputedObs,inputCol),
			from=0) 
	## extract x and y from the estimated density for the censored part of the sample
	estimatedDensityCensoredX <- 
			estimatedDensity$x[estimatedDensity$x < LOQ]
	estimatedDensityCensoredY <- 
			estimatedDensity$y[estimatedDensity$x < LOQ]
	## update the imputed observation by taking the 
	## conditional expectation
	updatedImputedObs <- sum(
			estimatedDensityCensoredX*estimatedDensityCensoredY)/
	sum(estimatedDensityCensoredY) 
	## Iterating updated imputed observation till it converge,
	## i.e., the difference between two 
	## computed conditional expectations becomes smaller
	## than epsilon
	# define number of iterations
	numIter <- 1
	while(abs(imputedObs-updatedImputedObs) > epsilon & numIter <= maxIter){ 
		imputedObs <- updatedImputedObs
		estimatedDensity <- density(c(imputedObs,inputCol), from=0)
		estimatedDensityCensoredX <- 
				estimatedDensity$x[estimatedDensity$x < LOQ]
		estimatedDensityCensoredY <- 
				estimatedDensity$y[estimatedDensity$x < LOQ]
		updatedImputedObs <- sum(
				estimatedDensityCensoredX*estimatedDensityCensoredY)/
		sum(estimatedDensityCensoredY)
		numIter <- numIter + 1
	}
	if (abs(imputedObs-updatedImputedObs) > epsilon){
		warning("Convergence is not achieved, increase maxIter or use a larger epsilon.")
	}
	return(updatedImputedObs)
}
### This function replaces the following function in the 
### original codes:
### kdens_iterations
