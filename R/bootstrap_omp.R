# void BLB_serial(double x[], double *result, double *gamma, int *s, int *R, int *n)
BLB_serial <- function(x, gamma, s = 15, r = 100) {
  ans <- .C("BLB_serial", as.double(x), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(x)))
  return(ans$result)
}


# void BLB_omp_on_s(double x[], double *result, double *gamma, int *s, int *R, int *n)
BLB_omp_on_s <- function(x, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_s", as.double(x), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(x)))
  return(ans$result)
}


# void BLB_omp_on_B(double x[], double *result, double *gamma, int *s, int *R, int *n)
BLB_omp_on_B <- function(x, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_B", as.double(x), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(x)))
  return(ans$result)
}