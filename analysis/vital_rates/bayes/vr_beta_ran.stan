
data { 

  // int<lower=1> nVital;      // N. of vital rates               
  int<lower=0> nYears;         // N. of years. Same for all vital rates.

  // Data for the growth model
  int<lower=0> nG;             // N. of data points for the growth model  
  int<lower=0> yearsG[nG];     // Index for years
  vector[nG] xG;               // log size at time t
  vector[nG] yG;               // log size at time t+1 
  
  // Data for the survival model
  int<lower=0> nS;             // N. of data points for the survival model  
  int<lower=0> yearsS[nS];     // Index for years
  vector[nS] xS;               // log size at time t
  vector[nS] xS2;              // log size2 at time t
  vector[nS] xS3;              // log size at time t
  int<lower=0,upper=1> yS[nS]; // Survival at time t+1. Values are either 0 or 1
  
  // Data for the flowering model
  int<lower=0> nF;             // N. of data points for the flowering model  
  int<lower=0> yearsF[nF];     // Index for years
  vector[nF] xF;               // log size at time t
  int<lower=0,upper=1> yF[nF]; // Flowering status at time t. Values are either 0or1
  
  // Data for the fertility model
  int<lower=0> nR;             // N. of data points for the fertility model  
  int<lower=0> yearsR[nR];     // Index for years
  vector[nR] xR;               // log size at time t
  int<lower=0> yR[nR];         // Number of flowers at time t

}

parameters {
  
  // means for intercept's random effects
  real yr_a_u_s;
  real yr_a_u_g;
  real yr_a_u_f;
  real yr_a_u_r;
  // means for slope's random effects
  real yr_b_u_s;
  real yr_b_u_g;
  real yr_b_u_f;
  real yr_b_u_r;
  
  // sd of intercept random effects
  real<lower=0> tau_a_s; 
  real<lower=0> tau_a_g; 
  real<lower=0> tau_a_f; 
  real<lower=0> tau_a_r; 
  // sd of slope random effects
  real<lower=0> tau_b_s; 
  real<lower=0> tau_b_g; 
  real<lower=0> tau_b_f; 
  real<lower=0> tau_b_r; 
  
  // placeholders of intercept random year effects
  real yr_a_s[nYears];
  real yr_a_g[nYears];
  real yr_a_f[nYears];
  real yr_a_r[nYears];

  // placeholders of slope random year effects
  real yr_b_s[nYears];
  real yr_b_g[nYears];
  real yr_b_f[nYears];
  real yr_b_r[nYears];

  real<lower=0> sigma_y;    // Residual st. dev. for growth model
  real b_s;               // Survival reg. slope
  real b_s2;              // Survival reg. slope
  real b_s3;              // Survival reg. slope
  real b_g;               // Growth reg. slope
  real b_f;               // Flowering reg. slope
  real b_r;               // Fertility reg. slope
  // real<lower=0> alphaR;  // Fertility dispersion parameter

}

model {

  // placeholders of indexes for the random year effects
  int indGY;   
  int indSY;
  int indFY;
  int indRY;
  // Placeholders for the surv(S), grow(G), flow(F) and fert(R) mods
  real mS[nS];
  real mG[nG]; 
  real mF[nF];
  real mR[nR];
  
  // Hyper-priors ()
  yr_b_u_s ~ normal(0, 20);
  yr_b_u_g ~ normal(0, 20);
  yr_b_u_f ~ normal(0, 20);
  yr_b_u_r ~ normal(0, 20);
  
  // hyperparameters
  for (n in 1:nYears){ // Random year effects 
    yr_a_s[n] ~ normal(yr_a_u_s, tau_a_s);
    yr_a_g[n] ~ normal(yr_a_u_g, tau_a_g);
    yr_a_f[n] ~ normal(yr_a_u_f, tau_a_f);
    yr_a_r[n] ~ normal(yr_a_u_r, tau_a_r);
    
    yr_b_s[n] ~ normal(yr_b_u_s, tau_b_s);
    yr_b_g[n] ~ normal(yr_b_u_g, tau_b_g);
    yr_b_f[n] ~ normal(yr_b_u_f, tau_b_f);
    yr_b_r[n] ~ normal(yr_b_u_r, tau_b_r);
  }

  // priors on parameters
  b_s  ~ normal(0, 100);   // Survival slope
  b_s2 ~ normal(0, 100);   // Survival slope
  b_s3 ~ normal(0, 100);   // Survival slope
  b_g  ~ normal(0, 100);    // Growth slope 
  b_f  ~ normal(0, 100);    // Flowering slope
  b_r  ~ normal(0, 100);    // Fertility slope 
  // alphaB ~ inv_gamma(0.001, 0.001);   // Fertility dispersion parameter 
  sigma_y ~ inv_gamma(0.001, 0.001);  // Growth model residual st. dev.   


  // Sampling ---------------------------------------------
  
  // survival
  for(nsurv in 1:nS){
    indSY     = yearsS[nsurv];
    mS[nsurv] = yr_a_s[indSY] + 
                b_s * xS[nsurv] +
                yr_b_s[indSY] * xS[nsurv] +
                b_s2 * xS2[nsurv] +
                b_s3 * xS3[nsurv];
  }
  yS ~ bernoulli_logit(mS);
  
  // growth
  for(ngrow in 1:nG){  // this (and subsequent 'for' loops) loops over random effects
    indGY     = yearsG[ngrow];
    mG[ngrow] = yr_a_g[indGY] + 
                b_g * xG[ngrow] +
                yr_b_g[indGY] * xG[ngrow];
  }  
  yG ~ normal(mG, sigma_y);

    
  // flowering
  for(nflow in 1:nF){
    indFY     = yearsF[nflow];
    mF[nflow] = yr_a_f[indFY] + 
                b_f * xF[nflow] +
                yr_b_f[indFY] * xF[nflow];
  }
  yF ~ bernoulli_logit(mF);
  
  // fertility
  for(nflow in 1:nR){
    indRY     = yearsR[nflow];
    mR[nflow] = yr_a_r[indRY] + 
                b_r * xR[nflow] + 
                yr_b_r[indRY] * xR[nflow];
  }
  yR ~ poisson_log(mR);

}
