data {
  int n;
  int<lower=0,upper=1> y[n];
  real x_size[n];
}

parameters {
  real b0;
  real b_size;
}

model {
  
  vector[n] m;
  
  for(ny in 1:n){
    m[ny] = b0 + b_size * x_size[ny];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha);
// }
