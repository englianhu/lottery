data {
  int N;
  int X[N, 3];
}

parameters {
  real<lower=0> lambda[6];
}

model {
  for(n in 1:N) {
      real lp[(1+min(X[n,1], X[n,2])) * (1+min(X[n,1], X[n,3])) * (1+min(X[n,2], X[n,3]))];
      int n_used;
      n_used = 0;
      for(y12 in 0:min(X[n,1], X[n,2]))
       for(y13 in 0:min(X[n,1], X[n,3]))
          for(y23 in 0:min(X[n,2], X[n,3]))    
            if((y12+y13) <= X[n,1] && (y12+y23) <= X[n,2] && (y13+y23) <= X[n,3]) {
              n_used = n_used + 1;
              lp[n_used] = poisson_lpdf(X[n,1]-y12-y13, lambda[1]) + 
                           poisson_lpdf(X[n,2]-y12-y23, lambda[2]) +
                           poisson_lpdf(X[n,3]-y13-y23, lambda[3]) +
                           poisson_lpdf(y12, lambda[4]) +
                           poisson_lpdf(y13, lambda[5]) +
                           poisson_lpdf(y23, lambda[6]);
            }
  
      increment_lpdf_prob(log_sum_exp(lp[1:n_used]));
  }
  
  lambda ~ gamma(0.1, 0.1);
}
