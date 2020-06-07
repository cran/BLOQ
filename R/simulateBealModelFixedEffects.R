#################################################################
### description:
###
### Simulating data from a Beal model with fixed effects only.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################



#' simulate data from Beal model with fixed effects
#' 
#' function to generate data from a Beal model with fixed effects
#' @encoding UTF-8
#' @param numSubjects scalar, number of subject which should be generated
#' @param clearance scalar, clearance
#' @param volumeOfDistribution scalar, volume of distribution
#' @param dose scalar, dose
#' @param timePoints vector of time points
#' @return generated sample with numSubjects as the number of rows
#' and length of timePoints as the number of columns
#' @details The model used to generate data at time t is as follows 
#' \deqn{y(t)=C(t)\exp(e(t)),}
#' where \eqn{C(t)}, the PK-model, is defined as follows:
#' \deqn{C(t) = \frac{\mathrm{dose}}{V_d} \exp{(CL.t)},}
#' with \eqn{V_d} the volume of distribution and \eqn{CL} as clearance.
#' The error model is consdiered as \eqn{e(t) \sim N(0, h(t))}, with:
#' \deqn{h(t) = 0.03 + 0.165  \frac{C(t)^{-1}}{C(1.5)^{-1} + C(t)^{-1}}} 
#' @seealso Beal S. L., Ways to fit a PK model with some data below 
#' the quantification limit, Journal of Pharmacokinetics
#' and Pharmacodynamics, 2001;28(\strong{5}):481–504.
#' @examples
#' set.seed(111)
#' simulateBealModelFixedEffects(10, 0.693,
#' + 		1, 1, seq(0.5,3,0.5))
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export
simulateBealModelFixedEffects <- function(numSubjects, clearance, volumeOfDistribution,dose, timePoints){
	simulateBealModelFixedEffects1Subject <- function(x, clearance, volumeOfDistribution,dose, timePoints){
		PKmodel<-function(timePoints){
			(dose/volumeOfDistribution)*exp(-clearance*timePoints)
		}
		errorVariance <- function(timePoints){
			0.03+0.165*((1/PKmodel(timePoints))/(1/(PKmodel(1.5))+(1/PKmodel(timePoints))))
		}
		errorModel <- function(timePoints){
			rnorm(1,0,(errorVariance(timePoints))^0.5)
		}
		computingBealModel <- function(timePoints){
			PKmodel(timePoints)*exp(errorModel(timePoints))
		}
		allTimePointsBealModel <- sapply(timePoints,computingBealModel)
		return(allTimePointsBealModel)
	}
	t(sapply(c(1:numSubjects), simulateBealModelFixedEffects1Subject, clearance, volumeOfDistribution, dose, timePoints))
}
