#################################################################
### description:
###
### Compute the log-censored likelihood of multivariate normal
### (additive error) or log-normal (multiplicative error)
### for a gievn dataset
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################


#' function to compute the cendored log-likelihood function
#' for a given dataset
#' @param parametersMVN vector parameters of the multivariate 
#' normal distribution
#' @param inputData matrix or data frame the input dataset
#' @param LOQ LOQ
#'@param isMultiplicative logical variable indicating whether
#' an addtive error model (FALSE- default) or a multiplicative
#' model should be used 
#' @return value of the censored multivariate normal 
#' log-likelihood
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @noRd
computeMVNLogLik <- function(parametersMVN, inputData, LOQ,
		isMultiplicative = FALSE){
	parametersLength <- length(parametersMVN)
	muLength <- (parametersLength+1)/3
	meanVec <- parametersMVN[1:(muLength)]
	covVec <- parametersMVN[(muLength+1):parametersLength]
	covMat <- convertCovandVec(covVec)
	invisible(ifelse(isMultiplicative,
					computedLogLikMVN <- 
							sum(apply(inputData,1,computeLogLikforEachRow, 
											meanVec, covMat, LOQ, 
											isMultiplicative = TRUE)),
					computedLogLikMVN <- 
							sum(apply(inputData,1,computeLogLikforEachRow, 
											meanVec, covMat, LOQ, 
											isMultiplicative = FALSE))))
	return(computedLogLikMVN)
}
### This function replaces the following functions in the 
### original codes:
### mvnormlik_cens_srs_off
### mvnormlik_cens_srs_off_g
