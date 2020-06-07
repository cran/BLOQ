#################################################################
### description:
###
### Implementing the imputation approach which simply replaces
### every BLOQ with a pre-defined constant. The two common
### options are using LOQ/2 or 0.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' imputing BLOQ's with a constant value
#' 
#' function to impute BLOQ observations by replacing them
#' with a constant value.
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar, limit of quantification value
#' @param constantValue scalar, the constant value which replaces
#' all BLOQ's, default is LOQ/2
#' @return the imputed dataset: a numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#' @examples
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#' + 		1, 1, seq(0.5,3,0.5))
#' # replacing BLOQ's with LOQ/2
#' imputeConstant(genDataFixedEffects, 0.1, 0.1/2)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
imputeConstant <- function(inputData, LOQ, constantValue){
	if (missing(constantValue)){
		constantValue <- LOQ/2
	}
	imputedData <- inputData
	imputedData[inputData<LOQ]=constantValue
	return(imputedData)
}
### This function replaces the following functions in the 
### original codes:
### replace_0
### replace_LOQ2