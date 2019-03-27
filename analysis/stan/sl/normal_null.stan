
data {
  int n_time;
  vector[n_time] y;
}

parameters {
  real alpha;
  real<lower=0> y_sd;
}

model {
  // model
  y ~ normal(alpha, y_sd);
}

generated quantities {
  vector[n_time] log_lik;
  
  for (n in 1:n_time)
    log_lik[n] = normal_lpdf(y[n] | alpha, y_sd);
}
