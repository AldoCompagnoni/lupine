
data { 

  // int<lower=1> nVital;      // N. of vital rates               
  int<lower=0> nYears;         // N. of years. Same for all vital rates.

  #Data for the growth model
  int<lower=0> nG;             // N. of data points for the growth model  
  int<lower=0> yearsG[nG];     // Index for years
  vector[nG] xG;               // log size at time t
  vector[nG] yG;               // log size at time t+1 
  
  #Data for the survival model
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
  
  // means for random effects
  real yr_u_s;
  real yr_u_g;
  real yr_u_f;
  real yr_u_r;
  
  // sd random effects
  real<lower=0> tau_s; 
  real<lower=0> tau_g; 
  real<lower=0> tau_f; 
  real<lower=0> tau_r; 
  
  # placeholders random year effects
  real yr_s[nYears];
  real yr_g[nYears];
  real yr_f[nYears];
  real yr_r[nYears];

  real<lower=0> sigma_y;    // Residual st. dev. for growth model
  real betaS;               // Survival reg. slope
  real betaS2;              // Survival reg. slope
  real betaS3;              // Survival reg. slope
  real betaG;               // Growth reg. slope
  real betaF;               // Flowering reg. slope
  real betaR;               // Fertility reg. slope
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
  
  // priors parameters
  betaS  ~ normal(0, 100);   #Survival slope
  betaS2 ~ normal(0, 100);   #Survival slope
  betaS3 ~ normal(0, 100);   #Survival slope
  betaG ~ normal(0, 100);   #Growth slope 
  betaF ~ normal(0, 100);   #Flowering slope
  betaR ~ normal(0, 100);   #Fertility slope 
  // alphaB ~ inv_gamma(0.001, 0.001);   #Fertility dispersion parameter 
  sigma_y ~ inv_gamma(0.001, 0.001);  #Growth model residual st. dev.   

  // hyperparameters
  for (n in 1:nYears){ #Random year effects 
    yr_s[n] ~ normal(yr_u_s, tau_s);
    yr_g[n] ~ normal(yr_u_g, tau_g);
    yr_f[n] ~ normal(yr_u_f, tau_f);
    yr_r[n] ~ normal(yr_u_r, tau_r);
  }

  // Sampling
  
  // survival
  for(nsurv in 1:nS){
    indSY <- yearsS[nsurv];
    mS[nsurv] <- yr_s[indSY] + 
                 betaS * xS[nsurv] +
                 betaS2 * xS2[nsurv] +
                 betaS3 * xS3[nsurv];
  }
  yS ~ bernoulli_logit(mS);
  
  // growth
  for(ngrow in 1:nG){  # this (and subsequent 'for' loops) loops over random effects
    indGY = yearsG[ngrow];
    mG[ngrow] = yr_g[indGY] + betaG * xG[ngrow];
  }  
  yG ~ normal(mG, sigma_y);

    
  // flowering
  for(nflow in 1:nF){
    indFY = yearsF[nflow];
    mF[nflow] = yr_f[indFY] + betaF * xF[nflow];
  }
  yF ~ bernoulli_logit(mF);
  
  // fertility
  for(nflow in 1:nR){
    indRY = yearsR[nflow];
    mR[nflow] = yr_r[indRY] + betaR * xR[nflow];
  }
  yR ~ poisson_log(mR);

}