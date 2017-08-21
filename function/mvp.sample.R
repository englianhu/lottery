#' mvp.sample
#'
#' @param theta.means underlying Poisson rates for the mean
#' @param theta.covs  underlying Poisson rates for the covariance
#' @param N           number of samples to be drawn
#'
#' @return            n-by-m matrix of n m-variate Poisson samples, where m is of size y1
#' @export
#'
#' @examples
mvp.sample <- function(theta.means, theta.covs, N = NULL) {
  if (is.null(N)) Y <- rbind(sample.matrix(theta.means), sample.vector(theta.covs, nrow(theta.means)))
  else Y <- sample.vector(c(theta.means, theta.covs), N)
  
  m <- if (is.null(N)) ncol(theta.means) else length(theta.means)
  A <- mvp.matrix(m)
  
  samples <- t(A %*% Y)
  return(samples)
}

sample.vector <- function(v, n) matrix(rpois(n * length(v), v), length(v), n)
sample.matrix <- function(m)    t(matrix(rpois(length(m), m), nrow(m), ncol(m)))

