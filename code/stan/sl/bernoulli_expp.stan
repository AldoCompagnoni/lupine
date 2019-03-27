
functions {
  real dexppow(real x, real mu, real sigma, real beta) {
    return((beta / (2 * sigma * tgamma(1.0/beta)) ) * exp(-(fabs(x - mu)/sigma)^beta));
  }
}

data {
  int n;
  int n_time;
  int n_lag;
  int n_site;
  int<lower=0> year_i[n];
  int<lower=0> site_i[n];
  int<lower=0,upper=1> y[n];
  matrix[n_time, n_lag] clim;
  real expp_beta;
}

parameters {
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=0,upper=n_lag> sens_sd;
  real beta_site[n_site];
  real beta;
}

transformed parameters {
  vector[n] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dexppow(i, sens_mu, sens_sd, expp_beta);
  
  sens = sens / sum(sens);
  
  for(i in 1:n)
    x[i] = sum(sens .* row(clim, year_i[i]));
}

model {
  
  real m[n];
  int s_i;   // placeholder for site identity 
  
  for(ny in 1:n){
    s_i    = site_i[ny];
    m[ny]  = beta_site[s_i] + beta * x[ny];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * x[year_i[ny]]);
// }
