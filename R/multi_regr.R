##Bootstrap function that performs regression of y on x and computes the stadard deviation of
# the regression coefficients

#' Bootstrap
#' @param data the original sample, as a vector
#' @param  data is a vector
#' @param B the number of bootstrap replications
#' @param n the number of resamples to generate
#' @return vector with the standard deviation of the regression coeffiecients
#' @examples
#' bootstrap(,,,,,)

 bootstrap <- function (x, y, B, n) {

  d <- ncol( x )
  b <- nrow( x )
  ans <- .C('bootstrap', as.double( x ), as.double( y ), as.double (numeric ( d )) , as.integer( b ), as.integer (B), as.integer( n ), as.integer( d ) )
  return (ans[[3]])

}

