data {
  int n;
  int n_year;
  int year_i[n];
  int<lower=0,upper=1> y[n];
}

parameters {
  real b_yr[n_year];    // intercept coeffs by year
}

model {
  
  int yr_i;
  real m[n];
  
  for(ny in 1:n){
    yr_i  = year_i[ny];
    m[ny] = b_yr[yr_i];
  }
  y ~ bernoulli_logit(m);

}

generated quantities {
  vector[n] log_lik;

  for (ny in 1:n)
    log_lik[ny] = bernoulli_logit_lpmf(y[ny] | b_yr[year_i[ny]]);
}
