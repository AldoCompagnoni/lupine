data {
  int n;
  int n_site;
  int site_i[n];
  int<lower=0,upper=1> y[n];
}

parameters {
  real beta[n_site];
}

model {
  
  real m[n];
  int s_i;   
  
  for(ny in 1:n){
    s_i = site_i[ny];
    m[ny] = beta[s_i];
  }
  y ~ bernoulli_logit(m);
}

// generated quantities {
//   vector[n] log_lik;
// 
//   for (ny in 1:n)
//     log_lik[ny] = bernoulli_logit_lpmf(y[ny] | beta[site_i[ny]]);
// }
