int_add <- function(a,b) {
  ans <- .C("int_add", as.integer(a), as.integer(b), res=as.integer(0))
  return(ans$res)
}

