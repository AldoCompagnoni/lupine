
data {
  int n;
  int n_time;
  int n_site;
  int<lower=0> year_i[n]; 
  int<lower=0> site_i[n];
  int<lower=0,upper=1> y[n];
  vector[n_time] clim_means;
  real expp_beta;
}

parameters {
  real beta;
  real beta_site[n_site];
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  int s_i;   // placeholder for site identity 
  
  for(ny in 1:n){
    y_i   = year_i[ny];
    s_i   = site_i[ny];
    m[ny] = beta_site[s_i] + beta * clim_means[y_i];
  }
  y ~ bernoulli_logit(m);
}

generated quantities {
  vector[n] log_lik;

  for (ny in 1:n)
    log_lik[ny] = bernoulli_logit_lpmf(y[ny] | beta_site[site_i[ny]] + beta * clim_means[year_i[ny]]);
}
