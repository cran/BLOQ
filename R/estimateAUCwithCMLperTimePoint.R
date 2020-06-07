#################################################################
### description:
###
### Implementing the approach which uses censored 
### maximum likelihood per time point to estimate time point
### specific mean and standard deviation which are needed
### to estimate AUC and its standard error.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' estimate AUC with censored maximum likelihood per time point
#'
#' function to estimate mean and standard error of each column
#' of data with BLOQ's using a censored maximum likelihood (CML) approach,
#' then use these estimates for estimating AUC and its standard
#' error
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
#' @return a list with three components: output of maxLik function,
#' estimated parameters for 
#' each column using censored maximum likelihood, and estimated
#' AUC and its standard error.
#' @seealso \href{https://www.rdocumentation.org/packages/maxLik/versions/1.3-4/topics/maxLik}{maxLik}
#' @examples
#' # generate data from Beal model with only fixed effects
#' set.seed(111)
#' genDataFixedEffects <- simulateBealModelFixedEffects(10, 0.693,
#'  		1, 1, seq(0.5,3,0.5))
#' # Multiplicative error model
#' estimateAUCwithCMLperTimePoint(genDataFixedEffects, 0.1, seq(0.5,3,0.5), TRUE)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
estimateAUCwithCMLperTimePoint <- function(inputData, LOQ, 
		timePoints, isMultiplicative = FALSE, onlyFitCML = FALSE,
		printCMLmessage = TRUE, optimizationMethod = NULL, CMLcontrol = NULL){
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
	## First based on the value of isMultiplicative the log of the
	## data as well as the LOQ will be taken.
	if (isMultiplicative){
		inputData <- log10(inputData)
		LOQ <- log10(LOQ)
	}
	## Cheking if there are any BLOQs, if not just jump to the
	## AUC calculations
	CML = NULL
	if (sum(inputData<LOQ, na.rm = TRUE) == 0){
		invisible(ifelse(isMultiplicative, 
				estCMLmean <- exp(apply(inputData, 2, mean)),
				estCMLmean <- apply(inputData, 2, mean)))
		estCMLsd <- apply(inputData, 2, sd)
	}else{
		#require ("maxLik")
		estCMLmean <- rep(0, length(timePoints))
		estCMLsd <- rep(0, length(timePoints))
		## Estimate the paramerers using censored maximum 
		## likelihood (CML) for each time point separately
		for (iCol in 1:ncol(inputData)){
			dataCol <- inputData[,iCol]
			# cheking missing values, as the censored maximum 
			# likelihood is estimated per time point, we ignore 
			# all missing value, but let the user know about it
			if (sum(is.na(dataCol))>0){
				if (sum(is.na(dataCol))==1){
					warning(paste("There is one missing value in column", 
									iCol, " that is ignored!"))
				}
				if (sum(is.na(dataCol))>1){
					warning(paste("There are", sum(is.na(dataCol)), 
									" missing values in column", iCol, 
									" that  are ignored!"))
				}
				idxObserved <- which(is.na(dataCol) == FALSE)
				dataCol <- dataCol[idxObserved]
			}
			censoredCol <- dataCol
			isCensored <- dataCol < LOQ
			if (any(isCensored)){
				censoredCol [dataCol < LOQ] <- LOQ
				# Defining the censored log-likelihood function
				logLikCensoredFun <- function(inputParameters){
					mean <- inputParameters[1]
					sd <- inputParameters[2]
					if(sd<0){
						return(NA)
					}
					sum(ifelse(isCensored, 
									pnorm(LOQ, mean = mean, sd = sd,
											log.p = TRUE), 
									dnorm(censoredCol, mean = mean, 
											sd = sd, log = TRUE)))
				}
				# estimate the parameters using CML
		if (is.null(optimizationMethod)){
			fitCML <- try(maxLik(logLik = logLikCensoredFun, 
							start = c(mean = 0, sd = 1), control = CMLcontrol))
		}else{
			fitCML <- try(maxLik(logLik = logLikCensoredFun, 
							start = c(mean = 0, sd = 1), method = optimizationMethod,control = CMLcontrol))
		}
				
				# controlling non-conbvergence issues: if the parameters
				# cannot be estimated using censored maximum likelihood
				# then the function stops
				CML[[iCol]] <- fitCML
				if (is.character(fitCML)){
						stop(paste("The optimization function stopped with the following message:",fitCML))
				}
					if (is.null(fitCML$estimate)){
					stop(paste("The function stopped because censored
maximum likelihood could not estimate the parameters for
 column", iCol))
				}
				if (printCMLmessage){
					print(paste("The message from maximizing censored log-likelihood", "for column", iCol,":", fitCML$message))
					flush.console()
				}	
				invisible(ifelse (isMultiplicative, estCMLmean [iCol] <- 
								exp(as.numeric(fitCML$estimate[1])), 
						estCMLmean [iCol] <- as.numeric(fitCML$estimate[1])))
				estCMLsd [iCol] <- as.numeric(fitCML$estimate[2])
			}else{
				estCMLmean[iCol] <- mean(dataCol)
				invisible(ifelse (isMultiplicative, 
								estCMLmean [iCol] <- exp(mean(dataCol)), 
								estCMLmean [iCol] <- mean(dataCol)))
				estCMLsd[iCol] <- sd(dataCol)
				CML[[iCol]] <- NULL
			}
		}
	}
	estCML <- cbind(estCMLmean,estCMLsd)
	colnames(estCML) <- c('mean','sd')
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
		
		estAUC <- sum(computedWeights * estCMLmean)
		invisible(ifelse(isMultiplicative,
						varAUC <- sum(((computedWeights*estCMLmean)^2)*
										(estCMLsd^2)),
						varAUC <- sum(((computedWeights*estCMLsd)^2))))
		estimatedAUCandVariance <- c(estAUC,sqrt(varAUC))
		names(estimatedAUCandVariance)=c('Estimated AUC','Std Err.')
	}else{
		estimatedAUCandVariance <- NULL
	}
	return(list(CML = CML, estCML=estCML,AUC=estimatedAUCandVariance))
}
### This function replaces the following functions in the 
### original codes:
### lik_data_means_gv
### lik_data_means_av