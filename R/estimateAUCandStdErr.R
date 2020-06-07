#################################################################
### description:
###
### Computing an estimate of area under the concentrations versus 
### time curve (AUC) and its variance in an PK-study with 
### non-compartmental (NCA) for additive amd multiplicative
### error models.
### approach.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################

#' Estimate AUC and its standard error
#' 
#' function to estimate AUC and compute standard error of this
#' estimate 
#' @encoding UTF-8
#' @param imputedData numeric matrix or data frame of size 
#' n by J (n the sample size and J the number of time points)
#' @param timePoints vector of time points
#' @param isMultiplicative logical variable indicating whether
#' an additive error model (FALSE) or a multiplicative error
#' model (TRUE) should be used 
#' @param na.rm logical variable indicating whether the rows with
#' missing values should be ignored or not. 
#' @return vector of length 2 with estimated AUC and its 
#' standard error
#' @examples 
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#' + 		1, 1, seq(0.5,3,0.5))
#' # Impute the data with BLOQ's with one of the provided methods,
#' # for example, here we use ROS
#' imputedDataROS <- imputeROS(genDataFixedEffects, 0.1)
#' # estimate AUC and its standard error
#' estimateAUCandStdErr(imputedDataROS,seq(0.5,3,0.5))
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
estimateAUCandStdErr <- function (imputedData,timePoints,
		isMultiplicative = FALSE, na.rm = FALSE){
	## The estimate of AUC and its variance are computed for
	## additive and multiplicative error models separately
	## With multiplicative error, dur to takign logarithm,
	## zero and negative values are not allowed, so first we check this 
	if (isMultiplicative & any(imputedData <= 0)){
		stop("Due to taking logarithm, negative or zero elements
 in the imputedData are not allowed with multiplicative 
error model.")
	}
	## The length of timePionts and number of columns of 
	## input Data should be the same.
	if (length(timePoints)!= ncol(imputedData)){
		stop("The length of timePionts and number of 
						columns of input Data should be the same.")
	}
	if (na.rm == FALSE & sum(is.na(imputedData)) > 0){
		stop("There are missing values in the data! You may sue na.rm = TRUE to remove rows with NA.")
	}
	if (na.rm == TRUE){
		imputedData <- imputedData[which(complete.cases(imputedData)),]
	}
	## First step is to compute the necessary weights
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
	
	## Having the weights, AUC estimate and its variance 
	## can be computed. 
	if (!isMultiplicative){
		meanCols <- colMeans(imputedData)
		estimatedAUC <- sum(computedWeights * meanCols)
		varEstimatedAUC <- var(imputedData%*%computedWeights)
	}
	if (isMultiplicative){
		meanCols <- colMeans(log(imputedData))
		estimatedAUC <- sum(computedWeights * exp(meanCols))
		varEstimatedAUC <- (computedWeights * exp(meanCols))%*%
				var(log(imputedData))%*%(computedWeights * exp(meanCols))
	}
	## Following general convention, standard error will be 
	## reported instead of the variance.
	estimatedAUCandStd <- c(estimatedAUC, 
			sqrt(varEstimatedAUC))
  names(estimatedAUCandStd)=c('Estimated AUC','Std Err.')
	return(estimatedAUCandStd)
}
### estimateAUCandStdErr replaces the following functrions 
### in the original codes:
### weig
### AUC_calc_single
### AUC_var_data
### AUC_val_var_g