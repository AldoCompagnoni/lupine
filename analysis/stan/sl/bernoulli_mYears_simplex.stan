
data {
  int n;
  int n_year;
  int K;
  int<lower=0> year_i[n]; 
  int<lower=0,upper=1> y[n];
  vector[n_year] clim1_means;
  vector[n_year] clim2_means;
  vector[n_year] clim3_means;
}

parameters {
  real b0;
  real b_c;
  simplex[K] theta_k;
}

transformed parameters {
  vector[n_year] x;
  
  for(i in 1:n_year){
    x[i] = (theta_k[1] * clim1_means[i]) + 
           (theta_k[2] * clim2_means[i]) + 
           (theta_k[3] * clim3_means[i]);
  }
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    y_i   = year_i[ny];
    m[ny] = b0 + b_c * x[y_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * clim_means[year_i[ny]]);
// }
