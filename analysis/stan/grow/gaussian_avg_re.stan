
data {
  int n;
  int n_year;
  int<lower=0> year_i[n]; 
  real y[n];
  real x_size[n];
  vector[n_year] clim_means;
}

parameters {
  real b0;
  real<lower=0> s_yr;
  
  real b_yr[n_year];
  real b_c;
  real b_size;
  real<lower=0> s;
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  // Hyperpriors
  b0    ~ normal(0,100);
  s_yr  ~ inv_gamma(0.001, 0.001);
  
  for(ii in 1:n_year){
    b_yr[ii] ~ normal(b0, s_yr);
  }

  for(ny in 1:n){
    y_i   = year_i[ny];
    m[ny] = b_yr[y_i] + b_size * x_size[ny] + b_c * clim_means[y_i];
  }
  y ~ normal(m,s);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * clim_means[year_i[ny]]);
// }
