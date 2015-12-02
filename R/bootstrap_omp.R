BLB_R <- function(data, gamma, FUN, ..., s=15, r=100) {
  n <- length(data)
  b <- round(n^gamma)
  xis <- numeric(s)
  
  subsamp_mat <- matrix(0, s, b)
  subsamp <- numeric(b)
  resample_mat <- matrix(0, r, n)
  
  for (i in 1:s) {
    subsamp <- sample(data, b)
    resample_mat <- matrix(sample(subsamp, r*n, replace = TRUE), r, n)
    theta <- apply(resample_mat, 1, FUN, ...)
    xis[i] <- sd(theta)
  }
  return(mean(xis))
}


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

# void BLB_omp_on_sB(double x[], double *result, double *gamma, int *s, int *R, int *n)
BLB_omp_on_sB <- function(x, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_sB", as.double(x), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(x)))
  return(ans$result)
}


# void bootstrap_for_clust(double x[], double *result, double *gamma, int *R, int *n)
bootstrap_for_clust_R <- function(i, x, gamma, r = 100) {
  # a is a dummy first argument, as can use with clusterApply
  ans <- .C("bootstrap_for_clust", as.double(x), result = as.double(0), as.double(gamma), as.integer(r), as.integer(length(x)), as.integer(i))
  return(ans$result)
}