
data {
  int n;
  int n_year;
  int M;
  int<lower=0> year_i[n];
  int<lower=0,upper=1> y[n];
  real x_size[n];
  real x_size2[n];
  matrix[M,n_year] clim1;
}

parameters {
  real b0;
  real<lower=0> s_yr;
  real b0_size;
  real<lower=0> s_size_yr;

  simplex[M] theta_m;
  real b_yr[n_year];
  real b_size_yr[n_year];
  real b_size2;
  real b_c;
}

transformed parameters {
  vector[n_year] x;
  
  for(i in 1:n_year)
    x[i] = sum(theta_m .* clim1[,year_i[i]]);
}

model {
  real m[n];
  int y_i;   // placeholder for year climate
  
  // Hyperpriors
  b0        ~ normal(0,100);
  s_yr      ~ inv_gamma(0.001, 0.001);
  b0_size   ~ normal(0,100);
  s_size_yr ~ inv_gamma(0.001, 0.001);
 
  for(ii in 1:n_year){
    b_yr[ii]      ~ normal(b0, s_yr);
    b_size_yr[ii] ~ normal(b0_size, s_size_yr);
  }

  for(ny in 1:n){
    y_i    = year_i[ny];
    m[ny]  = b_yr[y_i] + 
             b_size_yr[y_i] * x_size[ny] + 
             b_size2 * x_size2[ny] + 
             b_c * x[y_i];
  }
  y ~ bernoulli_logit(m);
}
