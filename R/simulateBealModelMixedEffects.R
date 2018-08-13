#################################################################
### description:
###
### Simulating data from a Beal model with fixed and random effects.
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
#' simulate data from Beal model with fixed and random effects
#' 
#' function to generate data from a Beal model with fixed effects
#' @encoding UTF-8
#' @param numSubjects scalar, number of subject which should be generated
#' @param clearance scalar, clearance
#' @param volumeOfDistribution scalar, volume of distribution
#' @param dose scalar, dose
#' @param varCompClearance scalar, standard error of
#' the normal distribution generating clearance
#' @param varCompVolumeOfDistribution scalar, standard error of
#' the normal distribution generating volume of distribution
#' @param timePoints vector of time points
#' @return generated sample with numSubjects as the number of rows
#' and length of timePoints as the number of columns
#' @details The model used to generate data at time t is as follows 
#' \deqn{y(t)=C(t)\exp(e(t)),}
#' where \eqn{C(t)}, the PK-model, is defined as follows:
#' \deqn{C(t) = \frac{\mathrm{dose}}{V_d} \exp{(CL.t)},}
#' with \eqn{V_d} the volume of distribution and \eqn{CL} as clearance.
#' The error model is consdiered as \eqn{e(t) \sim N(0, h(t))}, with:
#' \deqn{h(t) = 0.03 + 0.165  \frac{C(t)^{-1}}{C(1.5)^{-1} + C(t)^{-1}}.}
#' For the mixed effects model, \eqn{CL=\widetilde{CL} \exp{(\eta_1)}}, and
#' \eqn{V_d=\widetilde{V_d} \exp{(\eta_2)}}, where \eqn{\eta_1 \sim N(0, w_1^2)} and 
#' \eqn{\eta_1 \sim N(0, w_2^2)}. Note that \eqn{w_1} and \eqn{w_2} are specified by \emph{varCompClearance},
#'  and \emph{varCompVolumeOfDistribution} in the arguments, respectively.
#' @seealso Beal S. L., Ways to fit a PK model with some data below the quantification limit, Journal of Pharmacokinetics
#' and Pharmacodynamics, 2001;28(\strong{5}):481–504.
#' @examples
#' set.seed(111)
#' simulateBealModelMixedEffects(10, 0.693,
#' 		1, 1, 0.2,0.2, seq(0.5,3,0.5))
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @export

simulateBealModelMixedEffects <- function(numSubjects, clearance,
		volumeOfDistribution, dose, varCompClearance, 
		varCompVolumeOfDistribution, timePoints){
	simulateBealModelMixedEffects1Subject <- function(x,clearance,
			volumeOfDistribution,dose,varCompClearance, 
			varCompVolumeOfDistribution, timePoints){
		clearance <- clearance*exp(rnorm(1,0,varCompClearance))
		volumeOfDistribution<-volumeOfDistribution*exp(rnorm(1,0,
						varCompVolumeOfDistribution))
		PKmodel <- function(timePoints){
			(dose/volumeOfDistribution)*exp(-clearance*timePoints)
		}
		errorVariance <- function(timePoints){
			0.03+0.165*(
						(1/PKmodel(timePoints))/(1/(PKmodel(1.5))+
							(1/PKmodel(timePoints))))
		}
		errorModel <- function(timePoints){
			rnorm(1,0,(errorVariance(timePoints))^0.5)
		}
		computingBealModel <- function(timePoints){
			PKmodel(timePoints)*exp(errorModel(timePoints))
		}
		#browser()
		allTimePointsBealModel<-sapply(timePoints, 
				computingBealModel)
		return(allTimePointsBealModel)
	}
	t(sapply(c(1:numSubjects), 
					simulateBealModelMixedEffects1Subject, clearance, 
					volumeOfDistribution, dose, varCompClearance, 
					varCompVolumeOfDistribution, timePoints))
}
