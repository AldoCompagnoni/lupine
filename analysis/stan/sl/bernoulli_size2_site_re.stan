data {
  int n;
  int n_year;
  int n_site;
  int site_i[n];
  int<lower=0> year_i[n];
  int<lower=0> site_i[n];
  int<lower=0,upper=1> y[n];
  real x_size[n];
  real x_size2[n];
}

parameters {
  real b0;
  real<lower=0> s_yr;
  
  real b_yr[n_year];
  real b_size;
  real b_size2;
  real b_site[n_site];
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  int s_i;   // placeholder for site
  
  // Hyperpriors
  b0    ~ normal(0,100);
  s_yr  ~ inv_gamma(0.001, 0.001);
  
  for(ii in 1:n_year){
    b_yr[ii] ~ normal(b0, s_yr);
  }

  for(ny in 1:n){
    y_i   = year_i[ny];
    s_i   = site_i[ny];
    m[ny] = b_yr[y_i] + 
            b_site[s_i] +
            b_size * x_size[ny] +
            b_size2 * x_size2[ny];
  }
  y ~ bernoulli_logit(m);
  
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | beta[site_i[ny]]);
// }
