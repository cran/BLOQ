#################################################################
### description:
###
### Implementing the approach which uses full censored 
### maximum likelihood with multivariatre normal distribution 
### fitted in a pairwise fashion, to estimate time point specific 
### mean and standard deviation which are needed to estimate AUC 
### and its standard error. 
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' estimate AUCwith pairwise censored maximum likelihood
#' 
#' function to estimate mean and and covariance matrix of censored data using a 
#' full censored maximum likelihood approach via fitting all possible pairs, then use these 
#' estimates for estimating AUC and its standard error
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
#' censored maximum likelihood estimates will be calculated.
#' @param optimizationMethod single string specifying the method to be used for optimizing the log-likelihood, 
#' the default is NULL that allows the function to decide the about the best method. Otherwise, one can select among choices
#' available via R package maxLik: "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), 
#' "BFGSR" (for the BFGS algorithm implemented in R), 
#' "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), 
#' "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). 
#' Lower-case letters (such as "nr" for Newton-Raphson) are allowed.
#' @param CMLcontrol list of arguments to control 
#' convergence of maximization algorithm. It is the same argument
#' as control in the function maxLik in the R package maxLik
#' @param na.rm logical variable indicating whether the lines with missing values
#' should be ignored (TRUE, default) or not (FALSE). Note that, it will be applied 
#' for the sub-datasets regarding each pair.
#' @return a list with three components: output of maxLik function,
#' estimated parameters (mean vector and the covariance matrix) 
#' using censored maximum likelihood, and estimated
#' AUC and its standard error.
#' @seealso \href{https://www.rdocumentation.org/packages/maxLik/versions/1.3-4/topics/maxLik}{maxLik}
#' @examples
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#' 		1, 1, seq(0.5,3,1.5))
#' estimateAUCwithPairwiseCML(genDataFixedEffects, 0.1, seq(0.5,3,1.5))
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export

estimateAUCwithPairwiseCML <- function(inputData, LOQ, 
		timePoints, isMultiplicative = FALSE, onlyFitCML = FALSE,
		optimizationMethod = NULL,
		CMLcontrol = NULL, na.rm = TRUE){
	
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
	## missing data
	if (!na.rm & sum(apply(inputData, 2, is.na)) > 0){
		stop("When using censored maximum lielihood with
						multivariate normal, missing value are not allowed.")
	}
	
	# first define all possible pairs
	numTimePoints <- length(timePoints)
	allPairsIdx <- combn(numTimePoints, 2)
	numPairs <- ncol(allPairsIdx)
	estMu <- matrix(NA, numPairs, numTimePoints)
	estSigma <- matrix(NA, numPairs, numTimePoints^2)
	fitCML <- NULL
	for (iPairs in 1:numPairs){
		pairInputData <- inputData[,allPairsIdx[,iPairs]]
		# cheking the missing values and removing them if na.rm = TRUE
		if (na.rm){
			pairInputData <- pairInputData[complete.cases(pairInputData),]
		}
		# Here we check if for every subject, we have at least one 
		# non-censoted measurement. If both of the measurements are
		# censored, that particular subject is removed from computations
		# of the current pair, it will be announced by issuing a warnin message.
		numCensoredSubjectwise <- apply(pairInputData< LOQ, 1, sum)
		if (any(numCensoredSubjectwise >1)){
			pairInputData <- pairInputData[-which(numCensoredSubjectwise>1),]
			warning(paste("For pair ", paste(allPairsIdx[,iPairs],collapse="-")," the following subjects are removed due to being completely censored for this particular pair: ", which(numCensoredSubjectwise>1)))
		}
		fitPair <- estimateAUCwithFullCML(pairInputData, LOQ, 
				timePoints[allPairsIdx[,iPairs]], isMultiplicative, onlyFitCML = TRUE,
				printCMLmessage = FALSE, CMLcontrol, optimizationMethod = optimizationMethod,
				na.rm = FALSE)
		estMu[iPairs, allPairsIdx[,iPairs]] <- fitPair$estCML$mu
		tmpSigma <- matrix(NA, numTimePoints, numTimePoints)
		tmpSigma[allPairsIdx[,iPairs], allPairsIdx[,iPairs]] <- fitPair$estCML$Sigma
		estSigma [iPairs,] <- c(tmpSigma)
		fitCML[[iPairs]] <- fitPair$CML
	}
	pairwiseEstMu <- apply(estMu, 2, mean, na.rm = TRUE)
	pairwiseEstSigma <- matrix(apply(estSigma, 2, mean, na.rm = TRUE), numTimePoints, numTimePoints)
	estCML = list(mu = pairwiseEstMu, Sigma = pairwiseEstSigma)
	## Computing AUC and its variance
	if (!onlyFitCML){
		## Computing the necessary weights for calculating AUC
		## and its variance
		invisible(ifelse (numTimePoints >2,
						computedWeights <- c(timePoints[2]/2, 
								(timePoints[3:numTimePoints]
											-timePoints[1:(numTimePoints-2)])/2, 
								(timePoints[numTimePoints]-
											timePoints[(numTimePoints-1)])/2),
						computedWeights <- c(timePoints[2]/2, 
								(timePoints[numTimePoints]-
											timePoints[(numTimePoints-1)])/2)))
		estAUC <- sum(computedWeights * pairwiseEstMu)
		invisible(ifelse(isMultiplicative,
						varAUC <- (computedWeights * pairwiseEstMu) %*% 
								pairwiseEstSigma %*% (computedWeights * pairwiseEstMu) ,
						varAUC <- computedWeights %*% pairwiseEstSigma %*% 
								computedWeights))
		estimatedAUCandVariance <- c(estAUC,sqrt(varAUC))
		names(estimatedAUCandVariance)=c('Estimated AUC','Std Err.')
	}else{
		estimatedAUCandVariance <- NULL
	}
	return(list(pairwiseCML = fitCML, PairwiseEstCML=estCML,
					AUC=estimatedAUCandVariance))
}