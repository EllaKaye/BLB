
 bootstrap <- function (x, y, B, n) {

  d <- ncol( x )
  b <- nrow( x )
  ans <- .C('bootstrap', as.double( x ), as.double( y ), as.double (numeric ( d )) , as.integer( b ), as.integer (B), as.integer( n ), as.integer( d ) )
  return (ans[[3]])

}

