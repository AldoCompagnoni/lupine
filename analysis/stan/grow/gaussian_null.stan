
data {
  int n;
  real y[n];
  real x_size[n];
}

parameters {
  real b0;
  real b_size;
  real<lower=0> sigma;
}

model {
  real m[n];
  int y_i;   // placeholder for year climate
  
  for(ny in 1:n){
    m[ny]  = b0 + b_size * x_size[ny];
  }
  y ~ normal(m,sigma);
}
