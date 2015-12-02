#' Bootstraps for a regression coefficients
#' 
#' Bootstrap functions that performs regression of y on X and computes the standard deviation of the regression coefficients
#' 
#' The functions described here are designed as helper functions for the Bag of Little Bootstraps. When used in the BLB, either of these functions can be called on each of the s subsamples, with the matrix X and vector y being a subsample of the full data matrix and response vector.
#' 
#' \code{bootstrap_multidim} runs in serial over the r subsamples. 
#' @param X The data as a matrix
#' @param y The response variable for the regression
#' @param r the number of bootstrap resamples
#' @param n the number of data points in each resample (which is also the number of data points in the full data set)
#' @return vector with the standard deviation of the regression coeffiecients
#' @examples 
#' X <- matrix(rnorm(10), 5, 2)
#' y <- rnorm(5)
#' bootstrap_multidim(X, y, 1000, 30)
bootstrap_multidim <- function (X, y, r, n) {

  d <- ncol( X )
  b <- nrow( X )
  ans <- .C("bootstrap_multidim", as.double( X ), as.double( y ), result = as.double (numeric ( d )) , as.integer( b ), as.integer (r), as.integer( n ), as.integer( d ) )
  return(ans$result)
}

