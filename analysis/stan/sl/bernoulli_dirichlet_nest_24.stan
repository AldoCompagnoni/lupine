
data {
  int n;
  int n_year;
  int K;
  int M;
  int<lower=0> year_i[n];
  int<lower=0,upper=1> y[n];
  matrix[M,n_year] clim1;
  matrix[M,n_year] clim2;
}

parameters {
  simplex[K] theta_y;
  simplex[M] theta_m;
  real b0;
  real b_c;
}

transformed parameters {
  matrix[K,n_year] x_m;
  vector[n_year] x;
  
  for(i in 1:n_year){
    x_m[1,i] = sum(theta_m .* clim1[,year_i[i]]); 
    x_m[2,i] = sum(theta_m .* clim2[,year_i[i]]);
  }
  
  for(i in 1:n_year)
    x[i] = sum(theta_y .* x_m[,i]);
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    m[ny]  = b0 + b_c * x[y_i];
  }
  y ~ bernoulli_logit(m);
}
