
data {
  int n;
  int n_time;
  int n_lag;
  int<lower=0> year_i[n];
  int<lower=0,upper=1> y[n];
  matrix[n_lag, n_time] clim;
}

parameters {
  simplex[n_lag] theta;
  real alpha;
  real beta;
  //real beta_s[n_site];
  //real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  
  for(i in 1:n_time)
    x[i] = sum(theta .* clim[,year_i[i]]);
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    //s_i    = site_i[ny];
    m[ny]  = inv_logit(alpha + beta * x[y_i]);
  }
  y ~ bernoulli_logit(m);
}
