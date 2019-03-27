odata {
  int n;
  int n_year;
  int year_i[n];
  int<lower=0,upper=1> y[n];
}

parameters {
  
  real b0;              // mean intercept
  real b0_yr;              // mean intercept
  real<lower=0> sigma;  // intercept standard deviation

  real b_yr[n_year];    // intercept coeffs by year
  
}

model {
  
  real m[n];
  int yr_i;   
  
  // Hyperpriors
  b0    ~ normal(0,100);
  b0_yr ~ normal(0,100);
  sigma ~ inv_gamma(0.001, 0.001);
  
  for(ni in 1:n_year){ 
    b_yr[ni] ~ normal(b0_yr, sigma);
  }

  for(ny in 1:n){
    yr_i  = year_i[ny];
    m[ny] = b0 + b_yr[yr_i];
  }
  y ~ bernoulli_logit(m);

}

generated quantities {
  vector[n] log_lik;

  for (ny in 1:n)
    log_lik[ny] = bernoulli_logit_lpmf(y[ny] | b_yr[year_i[ny]]);
}
