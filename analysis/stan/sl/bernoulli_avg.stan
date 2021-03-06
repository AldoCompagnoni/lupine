
data {
  int n;
  int n_year;
  int<lower=0> year_i[n]; 
  int<lower=0,upper=1> y[n];
  vector[n_year] clim_means;
}

parameters {
  real b0;
  real b_c;
  // real beta_s[n_site];
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    y_i   = year_i[ny];
    m[ny] = b0 + b_c * clim_means[y_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * clim_means[year_i[ny]]);
// }
