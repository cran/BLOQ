#################################################################
### description:
###
### Implementing the imputation approach which uses Regression
### on order statistics (ROS) to impute the BLOQ's.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' imputing BLOQ's using regression on order statistics
#' 
#' function to impute BLOQ's with regression on order statistics
#' (ROS) approach.
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar limit of quantification value
#' @param isMultiplicative logical variable indicating whether
#' an additive error model (FALSE) or a multiplicative
#' model (TRUE) should be used 
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
#' imputeROS(genDataFixedEffects, 0.1)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
imputeROS <- function(inputData, LOQ, isMultiplicative = FALSE,
		useSeed = runif(1)){
	## Check whether the inputData includes any BLOQ
	## if not, the same dataset as inputData will be returned.
	imputedData <- inputData
	if (sum(imputedData<LOQ, na.rm = TRUE) == 0){
		print("There are not any BLOQ's in the dataset, so imputation is not needed")
	}else{
		for (iCol in 1:ncol(imputedData)){
			dataCol <- imputedData[,iCol]
			# Note that, as it could be possible to have
			# different number of measurement per time point
			# i.e., missing value, we ignore the NA's in the 
			# computation, and then replace them with NA's in the 
			# final imputed dataset
			idxObserved <- which(is.na(dataCol) == FALSE)
			dataCol <- dataCol[idxObserved]
			isCensored <- dataCol < LOQ
			if(any(isCensored)){
				invisible(ifelse(isMultiplicative,
						measuredObs <- sort(log(dataCol[!isCensored])),
						measuredObs <- sort(dataCol[!isCensored])))
				nObs <- length(dataCol)		
				nCensored <- sum((isCensored)) 
				nMeasured <- nObs-nCensored 
				if (nMeasured>1){
					## Calculations regarding regression on 
					## order statistics (ROS)
					empiricalExceedanceProbability <- 
							nMeasured/(nMeasured+nCensored)
					idxMeasured <- 1:nMeasured
					idxCensored <- 1:nCensored
					plottingPosiotionsForMeasuredObs <- 1-
							empiricalExceedanceProbability+
							empiricalExceedanceProbability*
							(idxMeasured/(nMeasured+1))
					plottingPosiotionsForCensoredObs <- 
							(idxCensored/(1+nCensored))*
							(1-empiricalExceedanceProbability)
					# Fitting the regression model of measured observations
					# on their corresponding plotting positions
					regressObsonQuantilesForMeasured <-
							lm (measuredObs~
											qnorm(plottingPosiotionsForMeasuredObs))
					imputingModelIntercept <- as.numeric(
							regressObsonQuantilesForMeasured$coefficients[1])
					imputingModelSlope <- as.numeric(
							regressObsonQuantilesForMeasured$coefficients[2])
					# Using the fitted model to compute imputed values
					# for censored observations (BLOQ's).
					tmpImputedCensoredObs <- imputingModelIntercept+
							imputingModelSlope*
							qnorm(plottingPosiotionsForCensoredObs)
					invisible(ifelse(isMultiplicative, imputedCensoredObs 
											<-	exp(tmpImputedCensoredObs), 
									imputedCensoredObs <- tmpImputedCensoredObs))
					# Ordering the imputed values according to the rules.
					if(iCol==1){
						# set seed to make the results reproducible
						set.seed(useSeed)
						imputingOrder <- sample(sum(isCensored))
					}else{
						imputingOrder <- rank(
								imputedData[, (iCol-1)][which(isCensored)])
					}
					dataCol[which(isCensored)] <- 
							sort(imputedCensoredObs)[imputingOrder]
				}else{
					print(paste('Due to only one measured observation, the imputation cannot be done for column ', iCol,' and the original values are kept'))
				}
				imputedData[idxObserved,iCol] <- dataCol
			}
		}
	}
	return(imputedData)
}
### This function replaces the following functions in the 
### original codes:
### ROS_impute_data_order
### ROS_impute_data_a_o