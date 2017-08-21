#' mvp.prob
#'
#' @param X     MVP sample
#' @param theta underlying Poisson rates
#' @param n.mc     number of Monte Carlo samples if method='MC"
#'
#' @return Estimate of the probability P(x) for a given vector theta of Poisson rates
#' @export
#'
#' @examples
mvp.prob <- function(X, theta, method, logarithm, n.mc = 10000) {
  switch(method,
         MC = {
           A <- mvp.matrix(length(X))
           probability <- sum(colSums((A %*% matrix(rpois(n.mc * length(theta), theta),
                                                    length(theta))) == X) == length(X)) / n.mc
           if (logarithm) probability <- log(probability)
         },
         recursive =  {
           if (is.null(dim(X))) {  # x is a single sample
             dict <<- array(-1, dim = X + 1)
             A <- mvp.matrix(length(X))
             
             if (logarithm) probability <- recur.prob.log(X, A, log(theta))
             else           probability <- recur.prob(X, A, theta)
           } else {               # x is an array of samples
             dict <<- array(-1, dim = apply(X, 2, max) + 1 + 1)
             A <- mvp.matrix(ncol(X))
             
             if (logarithm) probability <-  sum(apply(X, 1, function(x) recur.prob.log(X, A, log(theta))))
             else           probability <- prod(apply(X, 1, function(x) recur.prob(X, A, theta)))
           }
         },
         analytical = {
           A <- mvp.matrix(length(X))
           A2 <- A[, (length(X) + 1):ncol(A)]
           tmp <- X * A2 - (A2 - 1) * max(X);
           
           maximums <- apply(tmp, 2, min)
           y2 <- rep(0, length(maximums))
           probs <- rep(-Inf, prod(maximums + 1))
           
           j <- 1
           while (!(all(y2 == maximums))) {
             y1 <- X - A2 %*% y2
             if (sum(y1 < 0 ) == 0) probs[j] <- sum(dpois(c(y1, y2), theta, log = T))
             
             y2[1] <- y2[1] + 1
             for (i in 1:(length(y2) - 1)) {
               if (y2[i] > maximums[i]) {
                 y2[i] <- 0
                 y2[i + 1] <- y2[i + 1] + 1
               }
             }
             j <- j + 1
             
           }
           probs[j + 1] <- sum(dpois(c(X - A2 %*% y2, y2), theta, log = T))
           log.prob <- matrixStats::logSumExp(probs)
           
           if (logarithm) probability <- log.prob
           else           probability <- exp(log.prob)
         },
         stop('Unknown method. Please use "MC", "recursive", or "analytical"')
  )
  return(probability)
}


############################## PRIVATE HELPER FUNCTIONS ##############################

recur.prob.log <- function(X, A, logTheta) {
  if (sum(X < 0) > 0) return(-Inf)
  if (sum(X) == 0)    return(-exp(matrixStats::logSumExp(logTheta)))
  
  nnzf <- which.max(X > 0)
  indices <- which(A[nnzf, ] > 0)
  logp <- rep(-Inf, length(indices))
  for (i in 1:length(indices)) {
    idx <- indices[i]
    mvp <- matrix(X - A[, idx] + 1, 1)
    
    if (length(dict[mvp]) == 0) next
    if (dict[mvp] == -1) {
      result <- recur.prob.log(X - A[, idx], A, logTheta)
      dict[mvp] <<- result
    }
    logp[i] <- logTheta[idx] + dict[mvp]
  }
  return(matrixStats::logSumExp(logp) - log(X[nnzf]))
}

recur.prob <- function(X, A, theta) {
  if (sum(X < 0) > 0) return(0)
  if (sum(X) == 0) return(exp(-sum(theta)))
  
  nnzf <- which.max(X > 0)
  p <- 0
  for (idx in which(A[nnzf, ] > 0)) {
    mvp <- matrix(X - A[, idx] + 1, 1)
    
    if (length(dict[mvp]) == 0) next
    if (dict[mvp] == -1) {
      result <- recur.prob(X - A[, idx], A, theta)
      dict[mvp] <<- result
    }
    p <- p + theta[idx] * dict[mvp]
  }
  return(p / X[nnzf])
}

