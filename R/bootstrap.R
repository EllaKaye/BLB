samp_k_from_n <- function(k, n) {
  ans <- .C("samp_k_from_n", as.integer(k), as.integer(n), a = as.integer(numeric(k)))
  return(ans$a)
}

bootstrap <- function(x, B=1000) {
  ans <- .C("bootstrap", as.double(x), result = as.double(0), as.integer(B), as.integer(length(x)))
  return(ans$result)
}
