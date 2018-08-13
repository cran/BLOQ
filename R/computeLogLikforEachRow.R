#################################################################
### description:
###
### Compute the log-censored likelihood of multivariate normal
### (additive error) or log-normal (multiplicative error)
### for a given row of data
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' function to compute the cendored log-likelihood function
#' for a given row of the data
#' @param inputRow vector a given row of the data
#' @param meanVec vector mean 
#' @param covMat matrix covariance matrix
#' @param LOQ LOQ
#' @param isMultiplicative logical variable indicating whether
#' an addtive error model (FALSE- default) or a multiplicative
#' model should be used 
#' @return the value of log-likelihood function
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @noRd
computeLogLikforEachRow <- function(inputRow, meanVec, covMat, 
		LOQ, isMultiplicative = FALSE){
	isCensored <- inputRow < LOQ
	# require("mvtnorm")
	if(any(isCensored)){
		idxCensored <- which(isCensored %in% T)
		idxObserved <- which(isCensored %in% F)
		censoredCovMat <- covMat[idxCensored,idxCensored]
		observedCovMat <- covMat[idxObserved,idxObserved]
		censoredObservedCovMat <- 
				matrix(covMat[idxCensored,idxObserved],
						nrow=length(idxCensored))
		censoredMean <- meanVec[idxCensored]
		observedMean <- meanVec[idxObserved]
		conditionalMean <- as.vector(
				censoredMean + 
						censoredObservedCovMat%*%solve(observedCovMat)%*%
						(inputRow[idxObserved]-observedMean))
		conditionalCovMat <- censoredCovMat - 
				censoredObservedCovMat%*%
				solve(observedCovMat)%*%t(censoredObservedCovMat)
		# This defines the lower bound for normal (additive) or 
		# log-normal (multiplicative) errors
		invisible(ifelse(isMultiplicative, 
						lowerCensoredPart <- rep(-Inf,length(idxCensored)),
						lowerCensoredPart <- rep(0,length(idxCensored))))
		computedLogLik <- (log(as.numeric(
									pmvnorm(lower=lowerCensoredPart,upper=rep((LOQ)
					,length(idxCensored)),
			mean=conditionalMean,sigma=as.matrix(conditionalCovMat))))+
					log(dmvnorm(inputRow[idxObserved], 
									mean=observedMean,sigma=as.matrix(observedCovMat))))
	} else{
		computedLogLik<-log(dmvnorm(
						inputRow, mean=meanVec,sigma=as.matrix(covMat)))
	}
	return(computedLogLik)
}
### This function replaces the following functions in the 
### original codes:
### sumlog_function
### sumlog_function_g