# write model statement
cat(
  "
data {
  int<lower=1> n;                 // number of observations
  vector<lower=0>[n] Y[2];        // observations
  int length_u;                   // number of u indices to get inner summand
  vector<lower=0>[length_u] u;    // u indices
  int<lower=0> which_n[length_u]; // corresponding site indices
  matrix[n, length_u] u_mat;      // matrix of indicators to facilitate summation
}


parameters {
  vector<lower=0>[3] theta;       // expected values for component responses
}


transformed parameters {
  real theta_sum;
  vector[3] log_theta;
  vector[n] u_sum;
  vector[length_u] u_terms;


  theta_sum = sum(theta);
  log_theta = log(theta);


  for (i in 1:length_u){
    u_terms[i] = exp(binomial_coefficient_lpdf(Y[1, which_n[i]], u[i]) + 
                      binomial_coefficient_lpdf(Y[2, which_n[i]], u[i]) + 
                      lgamma(u[i] + 1) + 
                      u[i] * (log_theta[3] - log_theta[1] - log_theta[2]));
  }
  u_sum = u_mat * u_terms;
}


model {
  vector[n] loglik_el;
  real loglik;
  for (i in 1:n){
    loglik_el[i] = Y[1, i] * log_theta[1] - lgamma(Y[1, i] + 1) +  
                  Y[2, i] * log_theta[2] - lgamma(Y[2, i] + 1) + 
                  log(u_sum[i]);
  }
  loglik = sum(loglik_el) - n*theta_sum;
  increment_lpdf_prob(loglik);
}
  ", file="mvpois.stan")




# estimate parameters
library(rstan)
d = list(Y=Y, n=n, length_u=length_u, u=u, u_mat=u_mat, which_n=which_n)
init = stan("mvpois.stan", data = d, chains=0)
mod = stan(fit=init, 
            data=d, 
            iter=2000, chains=3, 
            pars=c("theta"))
traceplot(mod, "theta", inc_warmup=F)




# evaluate parameter recovery
post = extract(mod)
par(mfrow=c(1, 3))
for (i in 1:3){
  plot(density(post$theta[, i]))
  abline(v=theta[i], col="red")
}
par(mfrow=c(1, 1))


library(scales)
pairs(rbind(post$theta, theta), 
      cex=c(rep(1, nrow(post$theta)), 2), 
      pch=c(rep(1, nrow(post$theta)), 19), 
      col=c(rep(alpha("black", .2), nrow(post$theta)), "red"), 
      labels=c(expression(theta[1]), expression(theta[2]), expression(theta[3]))
      )
