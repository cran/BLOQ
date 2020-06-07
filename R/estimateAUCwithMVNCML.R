#################################################################
### description:
###
### Thsi function combines two methods of estimating AUC using 
### Multivatiate normal censored maximum likelihood: using the
### full likelihood with special structure for the covariance matruix
### or estimating it unstructured using a pairwise approach.
###
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' estimate AUC with multivariate normal censored maximum likelihood
#' 
#' function to estimate mean and and covariance matrix of censored data using a 
#' full censored maximum likelihood approach (with a 
#' special structure for the covariance matrix which only allows correlations between 
#' successive time points), then use these 
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
#' @param na.rm logical variable indicating whether the lines with missing values
#' should be ignored (TRUE, default) or not (FALSE).
#' @param isPairwise logical variable, if TRUE the unstructured covariance 
#' matrix will be estimated using pairwise approach, otherwise (FALSE, default)
#' the full maximum likelihood will be used with a special structure imposed on the covariance matrix.
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
#' estimateAUCwithMVNCML(genDataFixedEffects, 0.1, seq(0.5,3,1.5))
#' estimateAUCwithMVNCML(genDataFixedEffects, 0.1, seq(0.5,3,1.5),
#' isPairwise = TRUE)
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
estimateAUCwithMVNCML <- function(inputData, LOQ, 
		timePoints, isMultiplicative = FALSE, onlyFitCML = FALSE,
		printCMLmessage = TRUE, optimizationMethod = NULL,
		CMLcontrol = NULL, na.rm = TRUE, isPairwise = FALSE){
	## Use one of two options based on the value of isPairwise
	if (isPairwise){
		outputResults <- estimateAUCwithPairwiseCML (inputData, LOQ, 
				timePoints, isMultiplicative, onlyFitCML, optimizationMethod = optimizationMethod,
				CMLcontrol, na.rm)
	}else{
		outputResults <- estimateAUCwithFullCML (inputData, LOQ, 
				timePoints, isMultiplicative, onlyFitCML,
				printCMLmessage, optimizationMethod = optimizationMethod, CMLcontrol, na.rm)
	}
	return(outputResults)
}
