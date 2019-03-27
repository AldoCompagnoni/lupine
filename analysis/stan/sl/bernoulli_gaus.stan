
functions {
  real dnorm(real x, real mu, real sigma) {
    return((1 / sqrt(2 *pi()*pow(sigma, 2))) * exp(-((x-mu)^2) / (2*pow(sigma, 2))));
  }
}

data {
  int n;
  int n_time;
  int n_lag;
  int<lower=0> year_i[n];    
  int<lower=0,upper=1> y[n];
  matrix[n_time, n_lag] clim;
}

parameters {
  real<lower=0,upper=n_lag> sens_mu;
  real<lower=0,upper=50> sens_sd;
  real alpha;
  real beta;
  //real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dnorm(i, sens_mu, sens_sd);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    m[ny] = alpha + beta * x[y_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * x[year_i[ny]]);
// }

