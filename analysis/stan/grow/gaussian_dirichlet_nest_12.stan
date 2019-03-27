
data {
  int n;
  int n_year;
  int K;
  int M;
  int<lower=0> year_i[n];
  real y[n];
  real x_size[n];
  matrix[M,n_year] clim1;
}

parameters {
  simplex[K] theta_y;
  simplex[M] theta_m;
  real b0;
  real b_size;
  real b_c;
  real<lower=0> s;
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
    m[ny]  = b0 + b_size * x_size[ny] + b_c * x[y_i];
  }
  y ~ normal(m,s);
}
