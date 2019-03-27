
functions {
  real dweib(real x, real alphaw, real sigma) {
    return( (alphaw / sigma) * (x / sigma)^(alphaw-1) * exp(-((x/sigma)^alphaw)) );
  }
}

data {
  int n;
  int n_time;
  int n_site;
  int n_lag;
  int<lower=0> year_i[n];
  int<lower=0> site_i[n];
  int<lower=0,upper=1> y[n];
  matrix[n_time, n_lag] clim;
}

parameters {
  real<lower=0,upper=1000> alphaw;
  real<lower=0,upper=1000> sigma;
  real alpha;
  real beta;
  real beta_s[n_site];
  //real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dweib(i, alphaw, sigma);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  int s_i;   // placeholder for site identity
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    s_i    = site_i[ny];
    m[ny]  = alpha + beta * x[y_i] + beta_s[s_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * x[year_i[ny]]);
// }

