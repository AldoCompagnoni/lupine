
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
  
  // Data for the growth model
  int<lower=0> nG;             // N. of data points for the growth model  
  int<lower=0> yearsG[nG];     // Index for years
  int<lower=0> loc_g[nG];     // Index for location
  vector[nG] xG;               // log size at time t
  vector[nG] yG;               // log size at time t+1 
  
  // Data for the flowering model
  int<lower=0> nF;             // N. of data points for the flowering model  
  int<lower=0> yearsF[nF];     // Index for years
  int<lower=0> loc_f[nF];     // Index for location
  vector[nF] xF;               // log size at time t
  vector[nF] c_f;             // climate covariate
  int<lower=0,upper=1> yF[nF]; // Flowering status at time t. Values are either 0or1
  
  // Data for the fertility model
  int<lower=0> nR;             // N. of data points for the fertility model  
  int<lower=0> yearsR[nR];     // Index for years
  int<lower=0> loc_r[nR];     // Index for location
  vector[nR] xR;               // log size at time t
  vector[nR] c_r;             // climate covariate
  int<lower=0> yR[nR];         // Number of flowers at time t

}

parameters {
  
  // means for year intercept's random effects
  real yr_a_u_s;
  real yr_a_u_g;
  real yr_a_u_f;
  real yr_a_u_r;
  // means for year slope's random effects
  real yr_b_u_s;
  real yr_b_u_g;
  real yr_b_u_f;
  real yr_b_u_r;
  // means for location intercept's random effects
  real loc_a_u_s;
  real loc_a_u_g;
  real loc_a_u_f;
  real loc_a_u_r;
  // means for location slope's random effects
  real loc_b_u_s;
  real loc_b_u_g;
  real loc_b_u_f;
  real loc_b_u_r;
  
  // sd of year intercept random effects
  real<lower=0> yr_a_tau_s; 
  real<lower=0> yr_a_tau_g; 
  real<lower=0> yr_a_tau_f; 
  real<lower=0> yr_a_tau_r; 
  // sd of year slope random effects
  real<lower=0> yr_b_tau_s; 
  real<lower=0> yr_b_tau_g; 
  real<lower=0> yr_b_tau_f;
  real<lower=0> yr_b_tau_r; 
  // sd of location intercept random effects
  real<lower=0> loc_a_tau_s; 
  real<lower=0> loc_a_tau_g; 
  real<lower=0> loc_a_tau_f; 
  real<lower=0> loc_a_tau_r; 
  // sd of location slope random effects
  real<lower=0> loc_b_tau_s; 
  real<lower=0> loc_b_tau_g; 
  real<lower=0> loc_b_tau_f;
  real<lower=0> loc_b_tau_r; 
   
  // placeholders of year intercept random year effects
  real yr_a_s[nYears];
  real yr_a_g[nYears];
  real yr_a_f[nYears];
  real yr_a_r[nYears];
  // placeholders of year slope random year effects
  real yr_b_s[nYears];
  real yr_b_g[nYears];
  real yr_b_f[nYears];
  real yr_b_r[nYears];

  // placeholders of location intercept random year effects
  real loc_a_s[n_loc];
  real loc_a_g[n_loc];
  real loc_a_f[n_loc];
  real loc_a_r[n_loc];
  // placeholders of location slope random year effects
  real loc_b_s[n_loc];
  real loc_b_g[n_loc];
  real loc_b_f[n_loc];
  real loc_b_r[n_loc];

  real<lower=0> sigma_y;    // Residual st. dev. for growth model
  // real b_s;               // Survival reg. slope
  real b_s2;              // Survival reg. slope
  real b_s3;              // Survival reg. slope
  real b_c_s;             // Survival climate slope
  real b_c_f;             // Flowering climate slope
  real b_c_r;             // Fertility climate slope
  // real b_g;               // Growth reg. slope
  // real b_f;               // Flowering reg. slope
  // real b_r;               // Fertility reg. slope
  // real<lower=0> alphaR;  // Fertility dispersion parameter

}

model {

  // placeholders of indexes for the random YEAR effects
  int indGY;   
  int indSY;
  int indFY;
  int indRY;
  // placeholders of indexes for the random LOCATION effects
  int i_s_l;   
  int i_g_l;
  int i_f_l;
  int i_r_l;
  // Placeholders for the surv(S), grow(G), flow(F) and fert(R) mods
  real mS[nS];
  real mG[nG]; 
  real mF[nF];
  real mR[nR];
  
  // Hyper-priors (0 because slopes should be ~0)
  // yr_b_u_s ~ normal(0, 20);
  // yr_b_u_g ~ normal(0, 20);
  // yr_b_u_f ~ normal(0, 20);
  // yr_b_u_r ~ normal(0, 20);
  // loc_b_u_s ~ normal(0, 20);
  // loc_b_u_g ~ normal(0, 20);
  // loc_b_u_f ~ normal(0, 20);
  // loc_b_u_r ~ normal(0, 20);
  
  // hyperparameters years
  for (n in 1:nYears){ 
    yr_a_s[n] ~ normal(yr_a_u_s, yr_a_tau_s);
    yr_a_g[n] ~ normal(yr_a_u_g, yr_a_tau_g);
    yr_a_f[n] ~ normal(yr_a_u_f, yr_a_tau_f);
    yr_a_r[n] ~ normal(yr_a_u_r, yr_a_tau_r);
    
    yr_b_s[n] ~ normal(yr_b_u_s, yr_b_tau_s);
    yr_b_g[n] ~ normal(yr_b_u_g, yr_b_tau_g);
    yr_b_f[n] ~ normal(yr_b_u_f, yr_b_tau_f);
    yr_b_r[n] ~ normal(yr_b_u_r, yr_b_tau_r);
  }

  // hyperparameter locations
  for (n in 1:n_loc){ 
    loc_a_s[n] ~ normal(loc_a_u_s, loc_a_tau_s);
    loc_a_g[n] ~ normal(loc_a_u_g, loc_a_tau_g);
    loc_a_f[n] ~ normal(loc_a_u_f, loc_a_tau_f);
    loc_a_r[n] ~ normal(loc_a_u_r, loc_a_tau_r);
    
    loc_b_s[n] ~ normal(loc_b_u_s, loc_b_tau_s);
    loc_b_g[n] ~ normal(loc_b_u_g, loc_b_tau_g);
    loc_b_f[n] ~ normal(loc_b_u_f, loc_b_tau_f);
    loc_b_r[n] ~ normal(loc_b_u_r, loc_b_tau_r);
  }

  // priors on parameters
  // b_s  ~ normal(0, 100);   // Survival slope
  // b_s2 ~ normal(0, 100);   // Survival slope
  // b_s3 ~ normal(0, 100);   // Survival slope
  // b_g  ~ normal(0, 100);    // Growth slope 
  // b_f  ~ normal(0, 100);    // Flowering slope
  // b_r  ~ normal(0, 100);    // Fertility slope 
  // alphaB ~ inv_gamma(0.001, 0.001);   // Fertility dispersion parameter 
  sigma_y ~ inv_gamma(0.001, 0.001);  // Growth model residual st. dev.   


  // Sampling ---------------------------------------------
  
  // survival
  for(nsurv in 1:nS){
    indSY     = yearsS[nsurv];
    i_s_l     = loc_s[nsurv];
    mS[nsurv] = yr_a_s[indSY] + 
                yr_b_s[indSY] * xS[nsurv] +
                loc_a_s[i_s_l] + 
                loc_b_s[i_s_l] * xS[nsurv] +
                b_s2 * xS2[nsurv] +
                b_s3 * xS3[nsurv] + 
                b_c_s * c_s[nsurv];
  }
  yS ~ bernoulli_logit(mS);
  
  // growth
  for(ngrow in 1:nG){  // this (and subsequent 'for' loops) loops over random effects
    indGY     = yearsG[ngrow];
    i_g_l     = loc_g[ngrow];
    mG[ngrow] = yr_a_g[indGY]  + 
                yr_b_g[indGY]  * xG[ngrow] +
                loc_a_g[i_g_l] + 
                loc_b_g[i_g_l] * xG[ngrow];
  }  
  yG ~ normal(mG, sigma_y);

  // flowering
  for(nflow in 1:nF){
    indFY     = yearsF[nflow];
    i_f_l     = loc_f[nflow];
    mF[nflow] = yr_a_f[indFY] + 
                yr_b_f[indFY] * xF[nflow] +
                loc_a_f[i_f_l] + 
                loc_b_f[i_f_l] * xF[nflow] +
                b_c_f * c_f[nflow];
  }
  yF ~ bernoulli_logit(mF);
  
  // fertility
  for(nfert in 1:nR){
    indRY     = yearsR[nfert];
    i_r_l     = loc_r[nfert];
    mR[nfert] = yr_a_r[indRY]  + 
                yr_b_r[indRY]  * xR[nfert] +
                loc_a_r[i_r_l] + 
                loc_b_r[i_r_l] * xR[nfert] +
                b_c_r * c_r[nfert];
  }
  yR ~ poisson_log(mR);

}
