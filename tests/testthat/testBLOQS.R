#################################################################
### description:
###
### Testing the package BLOQ for several conditions
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################
context("Testing outcomes of functions in BLOQ!")


test_that ("posieiveAUCandSrdErr",{
			## This test checks whether the function which estimate AUC,
			## estimate it and its StdErr positive.
			set.seed(111)
			inputData <- simulateBealModelFixedEffects(10, 0.693,
					1, 1, seq(0.5,3,1.5))
			timePoints <- seq (0.5, 3, 1.5)
			AUCandStdErr <- estimateAUCandStdErr (inputData, timePoints)
			AUCandStdErrFullCML <- estimateAUCwithFullCML(inputData, 0.125, timePoints)$AUC
			AUCandStdErrCML <- estimateAUCwithCMLperTimePoint(inputData, 0.125, timePoints)$AUC
			expect_that(sum(AUCandStdErr<=0), equals(0))
			expect_that(sum(AUCandStdErrFullCML<=0), equals(0))
			expect_that(sum(AUCandStdErrCML<=0), equals(0))
			
		})


test_that ("nonImputedValues",{
			## This function checks if the imputing functions don't change
			## the non-BLOQ (shall we call them ALOQ ? :) above...) measurement. 
			set.seed(111)
			inputData <- simulateBealModelFixedEffects(10, 0.693,
					1, 1, seq(0.5,3,0.5))
			imputedCML <- imputeCML (inputData, 0.125)
			expect_that(sum(inputData[inputData >= 0.125] - imputedCML[inputData >= 0.125]), equals(0))
			imputedROS <- imputeROS (inputData, 0.125)
			expect_that(sum(inputData[inputData >= 0.125] - imputedROS[inputData >= 0.125]), equals(0))
			imputedConstant <- imputeConstant(inputData, 0.125, 0.125/2)
			expect_that(sum(inputData[inputData >= 0.125] - imputedConstant[inputData >= 0.125]), equals(0))
			imputedKernel <- imputeKernelDensityEstimation (inputData, 0.125, epsilon = 1e-05)
			expect_that(sum(inputData[inputData >= 0.125] - imputedKernel[inputData >= 0.125]), equals(0))
		})

test_that ("estimatedVariance",{
			## This function checks if the imputing functions don't change
			## the non-BLOQ (shall we call them ALOQ ? :) above...) measurement. 
			set.seed(111)
			inputData <- simulateBealModelFixedEffects(10, 0.693,
					1, 1, seq(0.5,3,1.5))
			timePoints <- seq (0.5, 3, 1.5)
			AUCandStdErrFullCML <- eigen(estimateAUCwithFullCML(inputData, 0.125, timePoints)$estCML$Sigma)$values
			AUCandStdErrCML <- estimateAUCwithCMLperTimePoint(inputData, 0.125, timePoints)$estCML[,2]
			expect_that(sum(AUCandStdErrCML<=0), equals(0))
			expect_that(sum(AUCandStdErrFullCML<=0), equals(0))
		})


