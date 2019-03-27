data {
  int n;
  int n_year;
  int n_site;
  int year_i[n_year];
  int site_i[n_year];
  int<lower=0,upper=1> y[n];
}

parameters {
  
  real b0_yr;            // mean intercept year
  real<lower=0> sig_yr;  // intercept standard deviation

  real b0_site;            // mean intercept site
  real<lower=0> sig_site;  // intercept standard deviation

  real b_yr[n_year];      // intercept coeffs by year
  real b_site[n_year];    // intercept coeffs by site
  
}

model {
  
  real m[n];
  int yr_i; 
  int st_i; 
  
  // Hyperpriors
  b0_yr   ~ normal(0,100);
  sig_yr  ~ inv_gamma(0.001, 0.001);
  
  b0_site ~ normal(0,100);
  sig_site~ inv_gamma(0.001, 0.001);
  
  // random effects
  for(ni in 1:n_year){ 
    b_yr[ni] ~ normal(b0_yr, sig_yr);
  }
  for(ns in 1:n_site){
    b_site[ns] ~ normal(b0_site, sig_site);
  }

  // full data
  for(ny in 1:n){
    yr_i  = year_i[ny];
    st_i  = year_i[ny];
    m[ny] = b_yr[yr_i] + b_site[st_i];
  }
  y ~ bernoulli_logit(m);

}

generated quantities {
  vector[n] log_lik;

  for (ny in 1:n)
    log_lik[ny] = bernoulli_logit_lpmf(y[ny] | b_yr[year_i[ny]]);
}
