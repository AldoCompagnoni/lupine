
data { 

  // int<lower=1> nVital;      // N. of vital rates               
  int<lower=0> N_YR;         // N. of years. Same for all vital rates.
  int<lower=0> N_LOC;         // N. of locations

  // Data for the survival model
  int<lower=0> N_S;             // N. of data points for the survival model  
  int<lower=0> yr_s[N_S];     // Index for years
  int<lower=0> loc_s[N_S];     // Index for location
  vector[N_S] xS;               // log size at time t
  vector[N_S] xS2;              // log size2 at time t
  vector[N_S] xS3;              // log size at time t
  vector[N_S] c_s;             // climate covariate
  int<lower=0,upper=1> yS[N_S]; // Survival at time t+1. Values are either 0 or 1
 
  // Data for the growth model
  int<lower=0> N_G;             // N. of data points for the growth model  
  int<lower=0> yr_g[N_G];     // Index for years
  int<lower=0> loc_g[N_G];     // Index for location
  vector[N_G] xG;               // log size at time t
  vector[N_G] yG;               // log size at time t+1 

  // Data for the flowering model
  int<lower=0> N_F;             // N. of data points for the flowering model  
  int<lower=0> yr_f[N_F];       // Index for years
  int<lower=0> loc_f[N_F];      // Index for location
  vector[N_F] xF;               // log size at time t
  vector[N_F] c_f;              // climate covariate
  int<lower=0,upper=1> yF[N_F]; // Flow status at time t. Values are either 0or1
 
  // Data for the fertility model
  int<lower=0> N_R;             // N. of data points for the fertility model  
  int<lower=0> yr_r[N_R];     // Index for years
  int<lower=0> loc_r[N_R];     // Index for location
  vector[N_R] xR;               // log size at time t
  vector[N_R] c_r;             // climate covariate
  int<lower=0> yR[N_R];         // Number of flowers at time t
 
}

parameters {
  
  // SURVIVAL ---------------------------------------
  
  // Hyper-parameters
  // means for year intercept's random effects
  real a_u_s;
  real b_u_s;
  // sd of year intercept random effects
  real<lower=0> a_tau_yr_s;
  real<lower=0> b_tau_yr_s; 
  real<lower=0> a_tau_loc_s;
  real<lower=0> b_tau_loc_s;

  // parameters
  vector[N_YR] a_zzz_yr_s;
  vector[N_YR] b_zzz_yr_s;
  vector[N_LOC] a_zzz_loc_s;
  vector[N_LOC] b_zzz_loc_s;
  
  real b_s2;    // Survival reg. slope
  real b_s3;    // Survival reg. slope
  real b_c_s;   // Survival climate slope


  // GROWTH ---------------------------------------
  
  // means for year intercept's random effects
  real a_u_g;
  real b_u_g;
  
  // Residual st. dev. for growth model
  real<lower=0> sigma_y;    
  
  // sd of year intercept random effects
  real<lower=0> a_tau_yr_g;
  real<lower=0> b_tau_yr_g; 
  real<lower=0> a_tau_loc_g;
  real<lower=0> b_tau_loc_g;

  // parameters
  vector[N_YR] a_zzz_yr_g;
  vector[N_YR] b_zzz_yr_g;
  vector[N_LOC] a_zzz_loc_g;
  vector[N_LOC] b_zzz_loc_g;
  
  // FLOWERING --------------------------------------
  
  // means for year intercept's random effects
  real a_u_f;
  real b_u_f;
  
  // sd of year intercept random effects
  real<lower=0> a_tau_yr_f;
  real<lower=0> b_tau_yr_f; 
  real<lower=0> a_tau_loc_f;
  real<lower=0> b_tau_loc_f;

  // parameters
  vector[N_YR] a_zzz_yr_f;
  vector[N_YR] b_zzz_yr_f;
  vector[N_LOC] a_zzz_loc_f;
  vector[N_LOC] b_zzz_loc_f;
  
  real b_c_f;   // Flowering climate slope
  
  
  // FERTILITY ---------------------------------------
  
  // Hyper-parameters
  // means for year intercept's random effects
  real a_u_r;
  real b_u_r;
  // sd of year intercept random effects
  real<lower=0> a_tau_yr_r;
  real<lower=0> b_tau_yr_r; 
  real<lower=0> a_tau_loc_r;
  real<lower=0> b_tau_loc_r;

  // parameters
  vector[N_YR] a_zzz_yr_r;
  vector[N_YR] b_zzz_yr_r;
  vector[N_LOC] a_zzz_loc_r;
  vector[N_LOC] b_zzz_loc_r;
  
  real b_c_r;   // Fertility climate slope
  
}

transformed parameters {
  
  // placeholders of year intercept random year effects
  vector[N_YR] a_yr_s;
  vector[N_YR] b_yr_s;
  vector[N_LOC] a_loc_s;
  vector[N_LOC] b_loc_s;
  
  // placeholders of year intercept random year effects
  vector[N_YR] a_yr_g;
  vector[N_YR] b_yr_g;
  vector[N_LOC] a_loc_g;
  vector[N_LOC] b_loc_g;

  // placeholders of year intercept random year effects
  vector[N_YR] a_yr_f;
  vector[N_YR] b_yr_f;
  vector[N_LOC] a_loc_f;
  vector[N_LOC] b_loc_f;

  // placeholders of year intercept random year effects
  vector[N_YR] a_yr_r;
  vector[N_YR] b_yr_r;
  vector[N_LOC] a_loc_r;
  vector[N_LOC] b_loc_r;

  // SURVIVAL ---------------------------------------
  a_yr_s  = a_u_s + a_tau_yr_s  * a_zzz_yr_s;
  b_yr_s  = b_u_s + b_tau_yr_s  * b_zzz_yr_s;
  a_loc_s = a_u_s + a_tau_loc_s * a_zzz_loc_s;
  b_loc_s = b_u_s + b_tau_loc_s * b_zzz_loc_s;

  // GROWTH -----------------------------------------
  a_yr_g  = a_u_g + a_tau_yr_g  * a_zzz_yr_g;
  b_yr_g  = b_u_g + b_tau_yr_g  * b_zzz_yr_g;
  a_loc_g = a_u_g + a_tau_loc_g * a_zzz_loc_g;
  b_loc_g = b_u_g + b_tau_loc_g * b_zzz_loc_g;

  // FLOWERING --------------------------------------
  a_yr_f  = a_u_f + a_tau_yr_f  * a_zzz_yr_f;
  b_yr_f  = b_u_f + b_tau_yr_f  * b_zzz_yr_f;
  a_loc_f = a_u_f + a_tau_loc_f * a_zzz_loc_f;
  b_loc_f = b_u_f + b_tau_loc_f * b_zzz_loc_f;

  // FERTILITY --------------------------------------
  a_yr_r  = a_u_r + a_tau_yr_r  * a_zzz_yr_r;
  b_yr_r  = b_u_r + b_tau_yr_r  * b_zzz_yr_r;
  a_loc_r = a_u_r + a_tau_loc_r * a_zzz_loc_r;
  b_loc_r = b_u_r + b_tau_loc_r * b_zzz_loc_r;

}

model {

  // Placeholders  
  
  // indexes for random effects
  int I_y_s;
  int I_l_s;
  int I_y_g;
  int I_l_g;
  int I_y_f;
  int I_l_f;
  int I_y_r;
  int I_l_r;

  // mean predictions
  real m_s[N_S];
  real m_g[N_G];
  real m_f[N_F];
  real m_r[N_R];

  // Survival ------------------------------------------------

  // Hyper-priors (0 because slopes should be close to it)
  a_u_s   ~ normal(0, 100);
  b_u_s   ~ normal(0, 100);

  a_tau_yr_s  ~ inv_gamma(0.001, 0.001);
  b_tau_yr_s  ~ inv_gamma(0.001, 0.001); 
  a_tau_loc_s ~ inv_gamma(0.001, 0.001); 
  b_tau_loc_s ~ inv_gamma(0.001, 0.001); 
  
  // random effects (year and location)
  a_zzz_yr_s  ~ normal(0, 10);
  b_zzz_yr_s  ~ normal(0, 10);
  a_zzz_loc_s ~ normal(0, 10);
  b_zzz_loc_s ~ normal(0, 10);
  
  // Growth ------------------------------------------------

  // Hyper-priors (0 because slopes should close to it)
  a_u_g   ~ normal(0, 100);
  b_u_g   ~ normal(0, 100);
  sigma_y ~ inv_gamma(0.001, 0.001);  // Growth model residual st. dev.   

  a_tau_yr_g  ~ inv_gamma(0.001, 0.001);
  b_tau_yr_g  ~ inv_gamma(0.001, 0.001); 
  a_tau_loc_g ~ inv_gamma(0.001, 0.001); 
  b_tau_loc_g ~ inv_gamma(0.001, 0.001); 
  
  // random effects (year and location)
  a_zzz_yr_g  ~ normal(0, 10);
  b_zzz_yr_g  ~ normal(0, 10);
  a_zzz_loc_g ~ normal(0, 10);
  b_zzz_loc_g ~ normal(0, 10);
  
  // Flowering ------------------------------------------------

  // Hyper-priors (0 because slopes should be close to it)
  a_u_f   ~ normal(0, 100);
  b_u_f   ~ normal(0, 100);

  a_tau_yr_f  ~ inv_gamma(0.001, 0.001);
  b_tau_yr_f  ~ inv_gamma(0.001, 0.001); 
  a_tau_loc_f ~ inv_gamma(0.001, 0.001); 
  b_tau_loc_f ~ inv_gamma(0.001, 0.001); 
  
  // random effects (year and location)
  a_zzz_yr_f  ~ normal(0, 10);
  b_zzz_yr_f  ~ normal(0, 10);
  a_zzz_loc_f ~ normal(0, 10);
  b_zzz_loc_f ~ normal(0, 10);
  
  // Fertility ------------------------------------------------

  // Hyper-priors (0 because slopes should be close to it)
  a_u_r   ~ normal(0, 100);
  b_u_r   ~ normal(0, 100);

  a_tau_yr_r  ~ inv_gamma(0.001, 0.001);
  b_tau_yr_r  ~ inv_gamma(0.001, 0.001); 
  a_tau_loc_r ~ inv_gamma(0.001, 0.001); 
  b_tau_loc_r ~ inv_gamma(0.001, 0.001); 
  
  // random effects (year and location)
  a_zzz_yr_r  ~ normal(0, 10);
  b_zzz_yr_r  ~ normal(0, 10);
  a_zzz_loc_r ~ normal(0, 10);
  b_zzz_loc_r ~ normal(0, 10);
  
  
  // Sampling ---------------------------------------------
  
  // survival
  for(II_s in 1:N_S){
    I_y_s     = yr_s[II_s];
    I_l_s     = loc_s[II_s];
    
    m_s[II_s] = a_yr_s[I_y_s]  + 
                b_yr_s[I_y_s]  * xS[II_s] +
                a_loc_s[I_l_s] +
                b_loc_s[I_l_s] * xS[II_s] +
                b_s2  * xS2[II_s] +
                b_s3  * xS3[II_s] + 
                b_c_s * c_s[II_s];
  }
  yS ~ bernoulli_logit(m_s);
  
  // growth
  for(II_g in 1:N_G){
    I_y_g     = yr_g[II_g];
    I_l_g     = loc_g[II_g];
    
    m_g[II_g] = a_yr_g[I_y_g]  + 
                b_yr_g[I_y_g]  * xG[II_g] +
                a_loc_g[I_l_g] +
                b_loc_g[I_l_g] * xG[II_g];
  }
  yG ~ normal(m_g, sigma_y);
  
  // flowering
  for(II_f in 1:N_F){
    I_y_f     = yr_f[II_f];
    I_l_f     = loc_f[II_f];
    
    m_f[II_f] = a_yr_f[I_y_f]  + 
                b_yr_f[I_y_f]  * xF[II_f] +
                a_loc_f[I_l_f] +
                b_loc_f[I_l_f] * xF[II_f] +
                b_c_f * c_f[II_f];
  }
  yF ~ bernoulli_logit(m_f);
  
  // flowering
  for(II_r in 1:N_R){
    I_y_r     = yr_r[II_r];
    I_l_r     = loc_r[II_r];
    
    m_r[II_r] = a_yr_r[I_y_r]  + 
                b_yr_r[I_y_r]  * xR[II_r] +
                a_loc_r[I_l_r] +
                b_loc_r[I_l_r] * xR[II_r] +
                b_c_r * c_r[II_r];
  }
  yR ~ poisson_log(m_r);
  
}
