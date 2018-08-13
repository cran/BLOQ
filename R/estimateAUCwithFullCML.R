#################################################################
### description:
###
### Implementing the approach which uses full censored 
### maximum likelihood with multivariatre normal distribution 
### (with a special strucutre for the covariance matrix) 
### to estimate time point specific mean and standard deviation 
### which are needed to estimate AUC and its standard error. 
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' estimate AUC with Full censored maximum likelihood
#' 
#' function to estimate mean and and covariance matrix of censored data using a 
#' full censored maximum likelihood approach (with a 
#' special structure for the covariance matrix which only allows correlations between 
#' successive time points), then use these 
#' estimates for estimating AUC and its standard error
#' 
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar, limit of quantification value
#' @param timePoints vector of time points
#' @param isMultiplicative logical variable indicating whether
#' an additive error model (FALSE) or a multiplicative error
#' model (TRUE) should be used 
#' @param onlyFitCML logical variable with FALSE as default, if TRUE only the
#' censored maximum likelihood estimates will be calculated
#' @param printCMLmessage logical variable with TRUE as default, if TRUE then
#' messages regarding the convergence status of censored 
#' log-likelihood maximization will be printed.
#' @param CMLcontrol list of arguments to control 
#' convergence of maximization algorithm. It is the same argument
#' as control in the function maxLik in the R package maxLik
#' @param na.rm logical variable indicating whether the lines with missing values
#' should be ignored (TRUE, default) or not (FALSE).
#' @return a list with three components: output of maxLik function,
#' estimated parameters (mean vector and the covariance matrix) 
#' using censored maximum likelihood, and estimated
#' AUC and its standard error.
#' @seealso \href{https://www.rdocumentation.org/packages/maxLik/versions/1.3-4/topics/maxLik}{maxLik}
#' @examples
#' #' # generate data from Beal model with only fixed effects
#' set.seed(123)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#'		1, 1, seq(0.5,3,1.5))
#' estimateAUCwithFullCML(genDataFixedEffects, 0.15, seq(0.5,3,1.5))
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
estimateAUCwithFullCML <- function(inputData, LOQ, 
		timePoints, isMultiplicative = FALSE, onlyFitCML = FALSE,
		printCMLmessage = TRUE, CMLcontrol = NULL, na.rm = TRUE){
	## With multiplicative error, dur to takign logarithm,
	## zero values are not allowed, so first we check this 
	if (isMultiplicative & any(inputData <= 0)){
		stop("Due to taking logarithm, negative or zero elements in 
the inputData are not allowed with multiplicative error model.")
	}
	## The length of timePionts and number of columns of 
	## input Data should be the same.
	if (length(timePoints)!= ncol(inputData)){
		stop("The length of timePionts and number of 
						columns of input Data should be the same.")
	}
	## For the multivariate normal missings are not allowed
	if (na.rm){
		inputData <- inputData[complete.cases(inputData),]
	}else{
		if (sum(apply(inputData, 2, is.na)) > 0){
			stop("When using censored maximum lielihood with
							multivariate normal, missing value are not allowed.")
		}
	}
	
	## First based on the value of isMultiplicative the log of the
	## data as well as the LOQ will be taken.
	if (isMultiplicative){
		inputData <- log(inputData)
		LOQ <- log(LOQ)
	}
	## computing starting values by first replacing all BLOQ's with
	## LOQ and then fit the likelihood to the NOW non-censored data
	startValuesData <- imputeConstant(inputData, LOQ, LOQ)
	#require("mvnmle")
	startValuesEst <- mlest(startValuesData)
	startMeanEst <- startValuesEst$muhat
	startCovMatEst <- startValuesEst$sigmahat
	startValues <- as.numeric(c(startMeanEst, 
					convertCovandVec(startCovMatEst)))
	## Computing the CMLE, the maximization method is fixed by
	## Conjugate Gradient
	#require("maxLik")
	invisible(ifelse(isMultiplicative,
					fitCML <- try(maxLik(logLik=computeMVNLogLik, 
							inputData = inputData, LOQ = LOQ, 
							isMultiplicative = TRUE, start = startValues, 
							control = CMLcontrol)),
					fitCML <- try(maxLik(logLik=computeMVNLogLik, 
							inputData = inputData, LOQ = LOQ, 
							isMultiplicative = FALSE, 
							start = startValues, 
							control = CMLcontrol))))
	
	if (is.character(fitCML)){
		stop(paste("The optimization function stopped with the following message:",fitCML))
	}
	if (is.null(fitCML$estimate)){
		stop("The function stopped because censored 
								maximum likelihood could not estimate the parameters!")
	}
	if (printCMLmessage){
		print(paste("The message from maximizing censored log-likelihood:", fitCML$message))
		flush.console()
	}		
	## extracting and converting estimated parameters
	muLength = length(timePoints)
	invisible(ifelse(isMultiplicative,
					estCMLmean <- exp(fitCML$estimate[1:muLength]),
					estCMLmean <- fitCML$estimate[1:muLength]))
	estCMLcovMat <- convertCovandVec(
			fitCML$estimate[(muLength+1):(3*muLength-1)])
	estCML <- list(mu = estCMLmean, Sigma = estCMLcovMat)
	## Computing AUC and its variance
	if (!onlyFitCML){
		## Computing the necessary weights for calculating AUC
		## and its variance
		numTimePoints <- length(timePoints)
		invisible(ifelse (numTimePoints >2,
				computedWeights <- c(timePoints[2]/2, 
						(timePoints[3:numTimePoints]
									-timePoints[1:(numTimePoints-2)])/2, 
						(timePoints[numTimePoints]-
									timePoints[(numTimePoints-1)])/2),
				computedWeights <- c(timePoints[2]/2, 
				(timePoints[numTimePoints]-
							timePoints[(numTimePoints-1)])/2)))
		estAUC <- sum(computedWeights * estCML$mu)
		invisible(ifelse(isMultiplicative,
						varAUC <- (computedWeights * estCML$mu) %*% 
								estCML$Sigma %*% (computedWeights * estCML$mu) ,
						varAUC <- computedWeights %*% estCML$Sigma %*% 
								computedWeights))
		estimatedAUCandVariance <- c(estAUC,sqrt(varAUC))
		names(estimatedAUCandVariance)=c('Estimated AUC','Std Err.')
	}else{
		estimatedAUCandVariance <- NULL
	}
	return(list(CML = fitCML, estCML=estCML,
					AUC=estimatedAUCandVariance))
}
### This function replaces the following functions in the 
### original codes:
### start_fr_off
### multi_norm_cens_srs_off
### multi_norm_cens_srs_off_g
### full_like_norm_off
### full_like_norm_mv_off
### full_like_log_off
### full_like_lognorm_mv_off


