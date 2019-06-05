
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
  
  // means for year intercept's random effects
  real yr_a_u_s;
   // means for year slope's random effects
  real yr_b_u_s;
  // means for location intercept's random effects
  real loc_a_u_s;
  // means for location slope's random effects
  real loc_b_u_s;

  // sd of year intercept random effects
  real<lower=0> yr_a_tau_s; 
  // sd of year slope random effects
  real<lower=0> yr_b_tau_s; 
  // sd of location intercept random effects
  real<lower=0> loc_a_tau_s; 
  // sd of location slope random effects
  real<lower=0> loc_b_tau_s; 

  // placeholders of year intercept random year effects
  real yr_a_s[nYears];
  // placeholders of year slope random year effects
  real yr_b_s[nYears];

  // placeholders of location intercept random year effects
  real loc_a_s[n_loc];
   // placeholders of location slope random year effects
  real loc_b_s[n_loc];

  real b_s2;              // Survival reg. slope
  real b_s3;              // Survival reg. slope
  real b_c_s;             // Survival climate slope

}

model {

  // placeholders of indexes for the random YEAR effects
  int indSY;   
  // placeholders of indexes for the random LOCATION effects
  int i_s_l;   
  // Placeholders for the surv(S), grow(G), flow(F) and fert(R) mods
  real mS[nS];

  // Hyper-priors (0 because slopes should be ~0)
  yr_a_u_s ~ normal(0, 100);
  yr_b_u_s ~ normal(0, 20);
  loc_a_u_s ~ normal(0, 100);
  loc_b_u_s ~ normal(0, 20);

  
  // hyperparameters years
  for (n in 1:nYears){ 
    yr_a_s[n] ~ normal(yr_a_u_s, yr_a_tau_s);
    yr_b_s[n] ~ normal(yr_b_u_s, yr_b_tau_s);
   }

  // hyperparameter locations
  for (n in 1:n_loc){ 
    loc_a_s[n] ~ normal(loc_a_u_s, loc_a_tau_s);
    loc_b_s[n] ~ normal(loc_b_u_s, loc_b_tau_s);
  }

  // priors on parameters
  b_s2 ~ normal(0, 100);   // Survival slope2
  b_s3 ~ normal(0, 100);   // Survival slope3
  b_c_s ~ normal(0, 100);   // Temperature slope
 
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
  
}
