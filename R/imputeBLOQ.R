#################################################################
### description:
###
### function to impute BLOQ's. The user can define column-specific
## methods to impute the BLOQ's.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################

#' impute BLOQ's with various methods
#' 
#' function to impute BLOQ's. The user can define column-specific methods to impute the BLOQ's.
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar, limit of quantification value
#' @param imputationMethod could be a single string or a vector of strings with the same length as 
#' the number of time points (ncol(inputData)). If it is left blank, then the imputation is done using 
#' kernel density estimation method for the columns with at least one non-BLOQ component. For all the
#' rest (only BLOQ) the constant imputation is used. The allowed values are 
#' "constant", "ros", "kernel", "cml" corresponding to constant imputation, 
#' imputing using regression on order statistics, imputing using kernel density estimator, and 
#' imputing using censored maximum likelihood, respectively. 
#' @param progressPrint logical variable indicating whether the imputation progress should be printed or not.
#' @param ... any other argument which should be changed according to the input arguments regarding 
#' the functions corresponding to different imputation methods.
#' @return a list with two components: imputed dataset, and the methods used to impute each column.
#' @examples
#' set.seed(111)
#' inputData <- simulateBealModelFixedEffects(10, 0.693,1, 1, seq(0.5,3,0.5))
#' LOQ = 0.125
#' imputeBLOQ(inputData, LOQ, 
#' 		imputationMethod = c("cml", "ros", "kernel","constant", "constant", "constant"), 
#' 		maxIter = 500, isMultiplicative = TRUE, constantValue = LOQ)
#' imputeBLOQ(inputData, LOQ, maxIter = 500, isMultiplicative = TRUE, 
#' constantValue = LOQ/5, epsilon = 1e-04)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export

imputeBLOQ <- function(inputData, LOQ, imputationMethod , progressPrint = FALSE, ...){
	## defining the method of imputation. If it's not specified, then the default
	## will be used. The default uses Kernel density estimator for columns
	## who have at least one non-BLOQ, and for only BLOQ columns imputeConstant will be used.
	## Otherwise, it can be either a single string or a vector of strings with the length
	## as the number of time points (number of columns in inputData). If former, then the same
	## imputation method is used for all columns, otherwise, each column will use the specified method.
	methodNames <- c("constant", "ros", "kernel", "cml")
	if (missing(imputationMethod)){
		imputationMethod <- rep("kernel", ncol(inputData))
		numBLOQ <- apply(inputData < LOQ, 2, sum, na.rm = TRUE)
		numObserved <- apply(apply(inputData, 2, complete.cases), 2, sum)
		imputationMethod[numBLOQ > numObserved-1] <- "constant"
	}else{
		if (sum(imputationMethod %in% methodNames) != length(imputationMethod)){
			stop("Please either set no imputationMethod 
or select it from 'contant', 'ros', 'kernel', 'cml'")
		}
		if (length(imputationMethod) == 1){
			imputationMethod <- rep(imputationMethod, ncol(inputData))
		}else{
			if (length(imputationMethod)!= ncol(inputData)){
				stop("If specified, the imputationMethod should be a single string, 
or a vector of strings with length the same as number of columns 
in the inputData.")
			}
		}
	}
	## making the input parameters for all functions
	inputParams <- list(...)
	inputParams[["LOQ"]] <- LOQ
	## defining set of input arguments for different methods
	paramsConstant <- c("inputData", "LOQ", "constantValue")
	paramsROS <- c("inputData", "LOQ", "isMultiplicative", "useSeed")
	paramsKernel <- c("inputData", "LOQ", "epsilon", "maxIter", "useSeed")
	paramsCML <- c("inputData", "LOQ", "isMultiplicative", "useSeed", "printCMLmessage", 
			"CMLcontrol")
	## Define an itnernal function to do the imputation based on selected method
	imputeMethod <- function(impMet){
		if (impMet == "constant"){
			inputParamsIDX <- which(names(inputParams) %in% paramsConstant)
			imputedData <- do.call(imputeConstant, inputParams[inputParamsIDX])
		}
		if (impMet == "ros"){
			inputParamsIDX <- which(names(inputParams) %in% paramsROS)
			imputedData <- do.call(imputeROS, inputParams[inputParamsIDX])
		}
		if (impMet == "kernel"){
			inputParamsIDX <- which(names(inputParams) %in% paramsKernel)
			imputedData <- do.call(imputeKernelDensityEstimation, inputParams[inputParamsIDX])
		}
		if (impMet == "cml"){
			inputParamsIDX <- which(names(inputParams) %in% paramsCML)
			imputedData <- do.call(imputeCML, inputParams[inputParamsIDX])
		}
		return(imputedData)
	}
	## If the same imputation method should be applied for all time points
	if (length(imputationMethod) == 1){
		inputParams[["inputData"]] <- inputData
		imputedData <- imputeMethod(imputationMethod)
	}else{
		imputedData <- matrix(0, nrow(inputData), ncol(inputData))
		for (iCol in 1:ncol(inputData)){
			if (progressPrint){
				print(paste("imputing column ", iCol, " with method ", imputationMethod[iCol]))
			}
			inputParams[["inputData"]] <- as.matrix(inputData[, iCol])
			imputedData[,iCol] <- imputeMethod(imputationMethod[iCol])
		}
	}
	return(list(imputedData = imputedData, imputationMethod = imputationMethod))
 }



 