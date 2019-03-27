
functions{
  real dgev(real y, real loc, real scale, real shape){
    
    real inv_scale;
    real neg_inv_shape; 
    real inv_shape_p1; 
    real x;
    real xx;
    real lp;

    x             = (y - loc)/scale;
    inv_scale     = 1/scale;
    neg_inv_shape = -1/shape;
    inv_shape_p1  = (1/shape) + 1;

    xx = 1 + shape * x;
    lp = log(inv_scale) - pow(xx,neg_inv_shape) - (inv_shape_p1 * log(xx));
     
    return exp(lp);
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
  real<lower=-2,upper=2> shape;
  real<lower=0,upper=n_lag> scale;
  real<lower=-(n_lag*2),upper=n_lag*2> loc;
  real alpha;
  real beta;
  //real beta_s[n_site];
  //real<lower=0> y_sd;
}

transformed parameters {
  vector[n_time] x;
  row_vector[n_lag] sens;
  
  for(i in 1:n_lag)
    sens[i] = dgev(i, loc, scale, shape);
  
  sens = sens / sum(sens);
  
  for(i in 1:n_time)
    x[i] = sum(sens .* row(clim, i));
}

model {
  
  real m[n];
  int y_i;   // placeholder for year climate
  //int s_i;   // placeholder for site identity
  
  for(ny in 1:n){
    y_i    = year_i[ny];
    //s_i    = site_i[ny];
    m[ny]  = alpha + beta * x[y_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | alpha + beta * x[year_i[ny]]);
// }

