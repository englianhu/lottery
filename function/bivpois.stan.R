## Simple bivariate poisson model in stan
## following parameterization in Karlis and Ntzoufras 2003
## Simulate data
n <- 50


# indpendent poisson components
theta <- c(2, 3, 1)


X_i <- array(dim=c(n, 3))
for (i in 1:3){
  X_i[,i] <- rpois(n, theta[i])
}


# summation to produce bivariate poisson RVs
Y_1 <- X_i[, 1] + X_i[, 3]
Y_2 <- X_i[, 2] + X_i[, 3]


Y <- matrix(c(Y_1, Y_2), 
            ncol=n, 
            byrow=T)


# Trick to generate a ragged array-esque set of indices for interior sum
# we want to use a matrix operation to calculate the inner sum 
# (from 0 to min(y_1[i], y_2[i])) for each of i observations
minima <- apply(Y, 2, min)
u <- NULL
which_n <- NULL
for (i in 1:n){
  indices <- 0:minima[i]
  u <- c(u, indices)
  which_n <- c(which_n, 
               rep(i, length(indices))
  )
}


length_u <- length(u) # should be sum(minima) + n


# construct matrix of indicators to ensure proper summation
u_mat <- array(dim=c(n, length_u))
for (i in 1:n){
  u_mat[i, ] <- ifelse(which_n == i, 1, 0)
}


# plot data
library(epade)
bar3d.ade(Y_1, Y_2, wall=3)

