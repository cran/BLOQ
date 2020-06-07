#################################################################
### description:
###
### Transformign a covariance matrix to a vecotr or vice versa
### Author: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
### maintainer: Vahid Nassiri <vahid.nassiri@openanalytics.eu>
#################################################################

#' function to covert a covariance matrix (vector) to a vector 
#' (covariance matrix) with a special structure allowing 
#' non-zero correlation only for two successive variables.
#' @param input matrix of vector which should be converted
#' @return converted matrix or vector
#' 
#' @author Vahid Nassiri, Helen Yvette Barnett
#' @noRd
convertCovandVec <- function (input){
	input <- as.matrix(input)
	## Check whether input is a matrix or a vector, so then we 
	## proceed to the appropriate command.
	if (1 %in% dim(input)){
		vecLength <- length(input)
		# This function solely cobnvert a covariance matrix, i.e.,
		# a positive definite symmetric matrix, into a vector within
		# the assumptions of this method, therefore, it is expected
		# that the lenfth of the given vector be an odd number
		if (vecLength /2 == floor(vecLength/2)){
			stop("The given vector should have an odd length!")
		}
		numCol <- 0.5*(vecLength+1)
		tmpOutput <- matrix(0,numCol,numCol)
		selectingMat <- row(tmpOutput) - col(tmpOutput)
		tmpOutput[(selectingMat <= 0) & (selectingMat >= -1)] <- 
				input
		# correcting negative diagonal elements
		tmpOutput[,which(diag(tmpOutput)<0)] <- 
				-1 * tmpOutput[,which(diag(tmpOutput)<0)] 
		# Get the original matrix from Cholesky decomposition
		output <- t(tmpOutput)%*%tmpOutput
	}else{
		# As only two successive time points are allowed to be
		# correlated, this is used to select such elements.
		selectingMat <- row(input) - col(input)
		output <- chol(input)[(selectingMat <= 0) & 
						(selectingMat >=-1)]
	}
	return(output)
}
### This function repalces the following functions in the
### original codes:
### cova2vec_off
### vec2cova_off
