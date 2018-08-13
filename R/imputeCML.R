#################################################################
### description:
###
### Implementing the imputation approach which uses censored 
### maximum likelihood per time point to estimate time point
### specific mean and standard deviation then use them to
### impute BLOQ's.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' imputing BLOQ's using censored maximum likelihood
#' 
#' function to impute BLOQ's using quantiles of a normal 
#' distribution with mean and standard error estimates using 
#' censored maximum likelihood
#' @encoding UTF-8
#' @param inputData numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#'  the input dataset
#' @param LOQ scalar, limit of quantification value
#' @param isMultiplicative logical variable indicating whether
#' an additive error model (FALSE) or a multiplicative error
#' model (TRUE) should be used 
#' @param useSeed scalar, set a seed to make the results 
#' reproducible, default is runif(1), it is used to randomly 
#' order the first imputed column (if the first column has any BLOQ's)
#' @param printCMLmessage logical variable with TRUE as default, if TRUE then
#' messages regarding the convergence status of censored 
#' log-likelihood maximization will be printed.
#' @param CMLcontrol list of arguments to control 
#' convergence of maximization algorithm. It is the same argument
#' as control in the function maxLik in the R package maxLik
#' @return the imputed dataset: a numeric matrix or data frame of the size 
#' n by J (n the sample size and J the number of time points)
#' @seealso \href{https://www.rdocumentation.org/packages/maxLik/versions/1.3-4/topics/maxLik}{maxLik}
#' @examples
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#' 		1, 1, seq(0.5,3,0.5))
#' imputeCML(genDataFixedEffects, 0.1, FALSE, 1)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
imputeCML <- function(inputData, LOQ, 
		isMultiplicative = FALSE, useSeed = runif(1), printCMLmessage = TRUE, 
		CMLcontrol = NULL){
	## Check whether the inputData includes any BLOQ
	## if not, the same dataset as inputData will be returned.
	imputedData <- inputData
	if (sum(imputedData<LOQ, na.rm = TRUE) == 0){
		print("There are not any BLOQ's in the dataset, so imputation is not needed")
	}else{
		## First we estimate the mean and stadard error of each
		## column using censored maximum likelihood
		estCML <- estimateAUCwithCMLperTimePoint(inputData, LOQ, 
				timePoints = rep(0,ncol(inputData)), isMultiplicative, onlyFitCML = TRUE,
				printCMLmessage, CMLcontrol)$estCML
		nCensored <- apply (inputData < LOQ, 2, sum)
		isCensored <- which(nCensored > 0)
		invisible(ifelse(isMultiplicative, 
						probBLOQ <- pnorm(log(LOQ),estCML[isCensored,1],
								estCML[isCensored,2]), 
						probBLOQ <- pnorm(LOQ,estCML[isCensored,1],
								estCML[isCensored,2])))
		for (iCensoredCol in 1:length(isCensored)){
			colCensored <- isCensored[iCensoredCol]
			## Distributing probability of being BLOQ among different
			## BLOQ positions and compute corresponding quantiles
			## using estimated mean and sd for each colums with
			## a BLOQ
			# find the quantiles as the imputed values, note that
			# as the estimate are found for multiplicative or additive
			# models via "estimateAUCwithCMLperTimePoint", the mean
			# and sd of the qnorm does not need to be transformed,
			# only at the end for multiplicative error one needs to
			# exponentiate the quantile.
			imputedValues <- qnorm(((1:nCensored[colCensored]) / 
								(nCensored[colCensored]+1)) * 
							probBLOQ[iCensoredCol], 
					mean = estCML[colCensored,1], 
					sd = estCML[colCensored,2])
			if (isMultiplicative){
				imputedValues = exp(imputedValues)
			}
			## Order the imputed values
			if(colCensored == 1){
				# set seed to make the results reproducible
				set.seed(useSeed)
				imputingOrder <- sample(nCensored[colCensored])
			}else{
				imputingOrder <- rank(imputedData[,(colCensored-1)]
								[which(inputData[,colCensored]< LOQ)])
			}
			imputedData[which(inputData[,colCensored]< LOQ), 
					colCensored]<- sort(imputedValues)[imputingOrder]
		}
	}
	return(imputedData)
}
### this function replaces the following function is the 
### original codes"
### impute_data3