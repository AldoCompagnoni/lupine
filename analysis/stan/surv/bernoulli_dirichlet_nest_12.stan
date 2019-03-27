
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
  simplex[M] theta_m;
  real b0;
  real b_size;
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
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    m[ny]  = b0 + b_size * x_size[ny] + b_size2 * x_size2[ny] + b_c * x[y_i];
  }
  y ~ bernoulli_logit(m);
}
