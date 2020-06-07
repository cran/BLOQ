#################################################################
### description:
###
### imputes the BLOQ observation using kernel density estimation.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################

#' imputing BLOQ's using kernel density estimation
#' 
#' function to impute BLOQ observations using kernel density 
#' estimation.
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar, limit of quantification value
#' @param epsilon scalar with 1e-07 as default, the difference between two iterations 
#' which achieving it would stop the procedure (convergence). 
#' @param maxIter scalar, the maximum number of iterations with 1000 as default.
#' @param useSeed scalar, set a seed to make the results 
#' reproducible, default is runif(1), it is used to randomly 
#' order the first imputed column (if the first column has any BLOQ's)
#' @return the imputed dataset: a numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#' @examples
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#' + 		1, 1, seq(0.5,3,0.5))
#' imputeKernelDensityEstimation(genDataFixedEffects, 0.1, epsilon = 1e-05)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
imputeKernelDensityEstimation <- function (inputData, LOQ, 
		epsilon = 1e-07, maxIter = 1000, useSeed = runif(1)){
	## Check whether the inputData includes any BLOQ
	## if not, the same dataset as inputData will be returned.
	imputedData <- inputData
	if (sum(imputedData<LOQ, na.rm = TRUE) == 0){
		print("There are not any BLOQ's in the dataset, so imputation is not needed")
	}else{
		for(iCol in 1:ncol(inputData)){
			observedCol <- inputData[,iCol] 
			imputedCol <- imputedData[,iCol] 
			idxObserved <- which(is.na(observedCol) == FALSE)
			dataCol <- observedCol[idxObserved]
			imputedCol <- imputedCol[idxObserved] 
			isCensored <- dataCol < LOQ
			if(any(isCensored)){
				numCensored <- sum(isCensored, na.rm = TRUE) 
				imputedValues <- rep(0,numCensored)
				## ordering the imputed values
				# if the very first column has any BLOQ then the order
				# is random
				if(iCol==1){ 
					set.seed(useSeed)
					# the following has been used as seed in the original codes
					#set.seed(numCensored)
					imputingOrder <- sample(1:numCensored)
				# for columns other than the first one the order is based 
				# on the previous column.
				}else{
					imputingOrder <- rank(imputedData[,
									(iCol-1)][which(isCensored)]) 
				}
				## imputing the first BLOQ
				imputedValues[1] <- 
						imputeKernelDensityEstimationInnerIteration(dataCol[!isCensored],
								LOQ, epsilon, maxIter) 
				## if there are more than one BLOQ, here they are imputed.
				if (numCensored > 1){
					for(iCensored in 2:numCensored){
						imputedValues[iCensored] <- 
									imputeKernelDensityEstimationInnerIteration(
										c(imputedValues[1:(iCensored-1)], dataCol[!isCensored]),
										LOQ, epsilon, maxIter) 
					}
				}
				imputedCol[which(isCensored)] <- 
						sort(imputedValues)[imputingOrder] 
				imputedData[idxObserved,iCol] <- imputedCol
			}
		}
	}
	return(imputedData)
}
### This function replaces the following function in the 
### original codes:
### kdens_iterations_multi_imp
