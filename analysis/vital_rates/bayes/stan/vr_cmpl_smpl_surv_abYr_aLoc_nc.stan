
data { 

  // int<lower=1> nVital;      // N. of vital rates               
  int<lower=0> nYears;         // N. of years. Same for all vital rates.
  int<lower=0> n_loc;         // N. of locations

  // Data for the survival model
  int<lower=0> nS;             // N. of data points for the survival model  
  int<lower=0> yearsS[nS];     // Index for years
  int<lower=0> loc_s[nS];     // Index for location
  vector[nS] xS;               // log size at time t
  vector[nS] xS2;              // log size2 at time t
  vector[nS] xS3;              // log size at time t
  vector[nS] c_s;             // climate covariate
  int<lower=0,upper=1> yS[nS]; // Survival at time t+1. Values are either 0 or 1
  
}

parameters {
  
  // Hyper-parameters
  
  // means for year intercept's random effects
  real a_u_s;
  real b_u_s;
  // sd of year intercept random effects
  real<lower=0> yr_a_tau_s;
  real<lower=0> yr_b_tau_s; 
  real<lower=0> loc_a_tau_s;
  real<lower=0> loc_b_tau_s;

  // parameters
  vector[nYears] yr_a_s_zzz;
  vector[nYears] yr_b_s_zzz;
  vector[n_loc]  loc_a_s_zzz;
  vector[n_loc]  loc_b_s_zzz;
  
  real b_s;     // Survival reg. slope
  real b_s2;    // Survival reg. slope
  real b_s3;    // Survival reg. slope
  real b_c_s;   // Survival climate slope

}

transformed parameters {
  
  // placeholders of year intercept random year effects
  vector[nYears] yr_a_s;
  vector[nYears] yr_b_s;
  vector[n_loc]  loc_a_s;
  vector[n_loc]  loc_b_s;
  
  yr_a_s  = a_u_s + yr_a_tau_s  * yr_a_s_zzz;
  yr_b_s  = b_u_s + yr_b_tau_s  * yr_b_s_zzz;
  loc_a_s = a_u_s + loc_a_tau_s * loc_a_s_zzz;
  loc_b_s = b_u_s + loc_b_tau_s * loc_b_s_zzz;
  
}

model {

  // Placeholders  
  
  // indexes for random effects
  int i_s_y;
  int i_s_l;
  // mean predictions
  real mS[nS];

  // Hyper-priors (0 because slopes should be ~0)
  a_u_s   ~ normal(0, 100);
  b_u_s   ~ normal(0, 100);

  yr_a_tau_s  ~ inv_gamma(0.001, 0.001);
  yr_b_tau_s  ~ inv_gamma(0.001, 0.001); 
  loc_a_tau_s ~ inv_gamma(0.001, 0.001); 
  loc_b_tau_s ~ inv_gamma(0.001, 0.001); 
  
  // random effects (year and location)
  yr_a_s_zzz  ~ normal(0, 10);
  yr_b_s_zzz  ~ normal(0, 10);
  loc_a_s_zzz ~ normal(0, 10);
  loc_b_s_zzz ~ normal(0, 10);
  
  // priors on single parameters
  b_s  ~ normal(0, 100);   // Survival slope2
  b_s2 ~ normal(0, 100);   // Survival slope2
  b_s3 ~ normal(0, 100);   // Survival slope3
  b_c_s ~ normal(0, 100);   // Temperature slope
 
  // Sampling ---------------------------------------------
  
  // survival
  for(nsurv in 1:nS){
    i_s_y     = yearsS[nsurv];
    i_s_l     = loc_s[nsurv];
    mS[nsurv] = yr_a_s[i_s_y] + 
                yr_b_s[i_s_y] * xS[nsurv] +
                loc_a_s[i_s_l] +
                loc_b_s[i_s_l] * xS[nsurv] +
                b_s2 * xS2[nsurv] +
                b_s3 * xS3[nsurv] + 
                b_c_s * c_s[nsurv];
  }
  yS ~ bernoulli_logit(mS);
  
}
