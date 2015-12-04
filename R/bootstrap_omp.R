#' Bag of Little Bootstraps for a vector
#' 
#' Implements Bag of Little Bootstraps on a 1-d data vector. The statistic which is bootstrapped is the mean, with the exception of BLB_R, which can accept a generic function as the statistic.
#' 
#' The functions described here implement the Bag of Little Bootstraps (BLB) with varying degrees of efficiency.
#' 
#' \code{BLB_R} implements BLB entirely in R. It is the most generic function, in that it can accept any function (which acts on a vector) as the statistic to be bootstrapped. It is also by far the slowest implementations.
#' 
#' The rest of the functions are all implemented in C, but are built in to only be able to bootstrap the mean.
#' 
#' \code{BLB_serial} does not take advantage of any parallelism. 
#' 
#' \code{BLB_omp_on_s} uses openMP to parallelise the for loop over the s subsampling steps, then runs a serial bootstrap function within each subsample.
#' 
#' \code{BLB_omp_on_r} does not parallelise over the s subsamples, but uses openMP to parallelise the for loop over the r resampling steps within each subsample.
#' 
#' \code{BLB_omp_on_rs} uses openMP to parallelise over both the subsampling and resampling loops.
#' 
#' \code{BLB_cluster} uses parallelism in R to send each of the s subsamples to different servers in the cluster, then uses openMP to parallelise the interations over the r resamples within each subsample.
#' 
#' @name BLB
#' @aliases BLB_serial
#' @aliases BLB_omp_on_s
#' @aliases BLB_omp_on_r
#' @aliases BLB_omp_on_sr
#' @aliases BLB_cluster
#' @param data A vector of the original sample
#' @param gamma Controls the size of the subsamples - each subsample is of size b = (floor(length(data) ^ gamma). Ideally in gamma is in [0.5, 1]
#' @param FUN A function to calculate the estimator of the parameter of interest
#' @param s The number of subsamples
#' @param r The number of bootstrap replications (resamples) per subsample
#' @return The mean across the s subsamples of the standard errors of the parameter estimates
#' @examples
#' X <- rnorm(1000000, 0, 5)
#' BLB_R(X, mean, gamma=0.6)
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

#' @rdname BLB
#' @examples 
#' BLB_serial(X, 0.9)
#' @export
BLB_serial <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_serial", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}
#' @rdname BLB
#' @examples 
#' BLB_omp_on_s(X, 0.6)
#' @export
BLB_omp_on_s <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_s", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

#' @rdname BLB
#' @examples 
#' BLB_omp_on_r(X, 0.6)
#' @export
BLB_omp_on_r <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_B", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

#' @rdname BLB
#' @examples 
#' BLB_omp_on_sr(X, 0.6)
#' @export
BLB_omp_on_sr <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_sB", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}


# void bootstrap_for_clust(double x[], double *result, double *gamma, int *R, int *n)
bootstrap_for_clust <- function(i, data, gamma, r = 100) {
  # a is a dummy first argument, as can use with clusterApply
  ans <- .C("bootstrap_for_clust", as.double(data), result = as.double(0), as.double(gamma), as.integer(r), as.integer(length(data)), as.integer(i))
  return(ans$result)
}

#' @rdname BLB
#' @examples 
#' cluster <- c("greywagtail", "greyheron", "greypartridge", "greyplover")
#' BLB_cluster(X, 0.6, cluster, s=16, r=96)
#' @export
BLB_cluster <- function(data, gamma, cluster, s=15, r=100) {
  clust <-makePSOCKcluster(cluster)
  clusterEvalQ(clust, library("BLB", lib.loc = "~/R"))
  lambda <- clusterApply(clust, 1:s, "bootstrap_for_clust", data=data, gamma=gamma, r=r)
  stopCluster(clust)  
  return(mean(unlist(lambda)))
}

# void BLB_md_omp_on_s(double x[], double y[], double *result, double *gamma, int *s, int *R, int *n, int *d)
#BLB_md_omp_on_s <- function(X, y, gamma, s=15, r=100) {
#  n <- nrow(X)
#  d <- ncol(X)
#  b <- floor(n^gamma)
#  ans <- .C("BLB_md_omp_on_s", as.double(t(X)), as.double(y), result = as.double(numeric(d)), as.double(gamma), as.integer(s), as.integer(r), as.integer(n), as.integer(d))
#  return(ans$result)
#}

# void BLB_serial_multinom(double x[], double *result, double *gamma, int *s, int *R, int *n)
BLB_serial_multinom <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_serial_multinom", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

BLB_omp_on_s_multinom <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_s_multinom", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

BLB_omp_on_r_multinom <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_B_multinom", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

BLB_omp_on_sr_multinom <- function(data, gamma, s = 15, r = 100) {
  ans <- .C("BLB_omp_on_sB_multinom", as.double(data), result = as.double(0), as.double(gamma), as.integer(s), as.integer(r), as.integer(length(data)))
  return(ans$result)
}

bootstrap_for_clust_multinom <- function(i, data, gamma, r = 100) {
  # a is a dummy first argument, as can use with clusterApply
  ans <- .C("bootstrap_for_clust_multinom", as.double(data), result = as.double(0), as.double(gamma), as.integer(r), as.integer(length(data)), as.integer(i))
  return(ans$result)
}


BLB_cluster_multinom <- function(data, gamma, cluster, s=15, r=100) {
  clust <-makePSOCKcluster(cluster)
  clusterEvalQ(clust, library("BLB", lib.loc = "~/R"))
  lambda <- clusterApply(clust, 1:s, "bootstrap_for_clust_multinom", data=data, gamma=gamma, r=r)
  stopCluster(clust)  
  return(mean(unlist(lambda)))
}




