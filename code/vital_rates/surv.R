setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(brms)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# set rstan options to parallel cores
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# data
lupine_08   <- read.csv("data/lupine_05_08.csv")
lupine_18   <- read.csv("data/lupine_08_18.csv")
clim        <- read.csv("data/lupine_fc_vars.csv")
enso        <- read.csv("data/enso_data.csv")
lupine_df   <- bind_rows(lupine_08,lupine_18)

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( stage_t0 != 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0) )

        
# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
m_obs     <- 7
m_back    <- 36
expp_beta <- 20 # this for the 

# format climate - need to select climate predictor first 
if(clim_var == 'prec'){
  clim_mat <- subset(clim, clim_var == "prate") %>%
                mutate( value = replace(value, value < 0, 0) ) %>% 
                month_clim_form("precip", years, m_back, m_obs)
} else{
  clim_var_input <- clim_var
  
  if( clim_var == 'airt') clim_input <- clim else clim_input <- enso
  
  clim_mat <- subset(clim_input, clim_var == clim_var_input) %>%
                month_clim_form(clim_var_input, years, m_back, m_obs)
}

# seedling data
surv_clim <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )
                # indices for STAN models
surv_df   <- surv_clim %>% 
                mutate( year_i = year %>% as.factor %>% as.numeric,
                        site_i = location %>% as.factor %>% as.numeric,
                        avgt0   = surv_clim %>% select(V1:V12) %>% rowSums,
                        avgtm1  = surv_clim %>% select(V13:V24) %>% rowSums,
                        avgtm2  = surv_clim %>% select(V25:V36) %>% rowSums )

# climate data
clim_pred <- dplyr::select(surv_df, year) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select(-year)


# fit GLMM models --------------------------------------------------------

mod1 <- glmer(surv_t1 ~ log_area_t0 * avgt0 + (1 | year), data= surv_df, 
              family = binomial )
mod2 <- glmer(surv_t1 ~ log_area_t0 + avgtm1 + (1 | year), data= surv_df, 
              family = binomial )
mod3 <- glmer(surv_t1 ~ log_area_t0 + avgtm2 + (1 | year), data= surv_df, 
              family = binomial )

# write it out
bind_rows( fixef(mod1),fixef(mod2), fixef(mod3) ) %>% 
  write.csv(paste0('results/lme4/','surv_',clim_var,'.csv'),row.names=F)


# mod1 <- glmer(surv_t1 ~ log_area_t0 * avgt0  + (1 | year), data= surv_df, 
#               family = binomial )
# mod2 <- glmer(surv_t1 ~ avgt0 + (log_area_t0 | year), data= surv_df, 
#               family = binomial )
# mod3 <- glmer(surv_t1 ~ avgt0 + (1 | year) + (log_area_t0|location), data= surv_df, 
#               family = binomial )
# mod4 <- glmer(surv_t1 ~ avgt0 + (log_area_t0 | year) + (log_area_t0|location), data= surv_df, 
#               family = binomial )
# AIC(mod1,mod2,mod3,mod4)


# fit stan models ----------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n      = nrow(surv_df),
  n_year = surv_df$year_i %>% unique %>% length,
  yr_bck = m_back / 12,
  n_site = surv_df$site_i %>% unique %>% length,
  n_lag  = ncol(clim_pred),
  y      = surv_df$surv_t1,
  x_size = surv_df$log_area_t0,
  clim   = clim_pred,
  clim_means = rowMeans(clim_pred),
  year_i = surv_df$year_i,
  site_i = surv_df$site_i,
  expp_beta = expp_beta,
  
  # climate variables
  clim1  = t(clim_pred)[1:12 ,],
  clim2  = t(clim_pred)[13:24,],
  clim3  = t(clim_pred)[25:36,],
  clim1_means = rowMeans( clim_pred[,1:12] ),
  clim2_means = rowMeans( clim_pred[,13:24] ),
  clim3_means = rowMeans( clim_pred[,25:36] ),
  K      = ncol(clim_pred) / 12,
  M      = 12
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter   = 4000, 
  thin   = 2, 
  chains = 4
)

# # Average of previous 3 years
# fit_avg <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_avg.stan"),
#   data = dat_stan,
#   pars = c('b0', 'b_size', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # store results
# out_summ  <- summary(fit_avg)$summary
# out_post  <- extract(fit_avg) %>% as.data.frame
# 
# write.csv(out_summ, paste0('results/surv_a_summ_',clim_var,'.csv'))
# write.csv(out_post, paste0("results/surv_a_post_",clim_var,".csv"),row.names=F)
# 
# 
# # Multiple years, weighted with simplex
# fit_3yr <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_mYears_simplex.stan"),
#   data = dat_stan,
#   pars = c('theta_k', 'b0', 'b_size', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # store results
# out_summ  <- summary(fit_3yr)$summary
# out_post  <- extract(fit_3yr) %>% as.data.frame
# 
# write.csv(out_summ, paste0("results/surv_sam_3y_summ_",clim_var,".csv"))
# write.csv(out_post, paste0("results/surv_sam_3y_post_",clim_var,".csv"),row.names=F)


# Average of previous 3 years
fit_avg_re <- stan(
  file = paste0("analysis/stan/surv/bernoulli_avg_re.stan"),
  data = dat_stan,
  pars = c('b0', 's_yr', 'b_yr', 'b_size', 'b_c'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# store results
out_summ  <- summary(fit_avg_re)$summary
out_post  <- extract(fit_avg_re) %>% as.data.frame

write.csv(out_summ, paste0('results/surv_a_summ_re_',clim_var,'.csv'))
write.csv(out_post, paste0("results/surv_a_post_re_",clim_var,".csv"),row.names=F)


# Multiple years, weighted with simplex
fit_3yr_re <- stan(
  file = paste0("analysis/stan/surv/bernoulli_mYears_simplex_re.stan"),
  data = dat_stan,
  pars = c('theta_k', 'b0', 's_yr', 'b_yr', 'b_size', 'b_c'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# store results
out_summ  <- summary(fit_3yr_re)$summary
out_post  <- extract(fit_3yr_re) %>% as.data.frame

write.csv(out_summ, paste0("results/surv_sam_3y_summ_re_",clim_var,".csv"))
write.csv(out_post, paste0("results/surv_sam_3y_post_re_",clim_var,".csv"),row.names=F)


# # 12 month dirichlet
# fit_12_nest <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest_12.stan"),
#   data = dat_stan,
#   pars = c('theta_m', 'b0', 'b_size', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.8, stepsize = 0.001, max_treedepth = 10)
# )
# 
# # store results
# out_summ  <- summary(fit_12_nest)$summary
# out_post  <- extract(fit_12_nest) %>% as.data.frame
# 
# write.csv(out_summ, paste0("results/surv_sam_12m_summ_",clim_var,".csv"))
# write.csv(out_post, paste0("results/surv_sam_12m_post_",clim_var,".csv"),row.names=F)


# fit1_m  <- brm(area_t1 ~ area_t0, data = surv_df )
# fit1    <- brm(log_area_t1 ~ log_area_t0 + (1|year), 
#                data = surv_df %>% mutate(year=as.factor(year)) )
# fit_f   <- brm(surv_t1 ~ year, 
#                data = mutate(sl_df, year = as.factor(year) ), 
#                family = bernoulli(link = "logit") )
# 
# 
# # null model
# fit_mean <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_null.stan"),
#   data = dat_stan,
#   pars = c('b0', 'b_size'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # random year effect only
# fit_year <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_ran_yr.stan"),
#   data = dat_stan,
#   pars = c('b0', 's_yr', 'b_yr', 'b_size'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # update data list
# dat_stan$clim1        <- t(clim_pred)[1:12 ,]
# dat_stan$clim2        <- t(clim_pred)[13:24,]
# dat_stan$clim3        <- t(clim_pred)[25:36,]
# dat_stan$K            <- ncol(dat_stan$clim) / 12
# dat_stan$M            <- 12
# 
# # power exponential moving window
# fit_36_nest <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest.stan"),
#   data = dat_stan,
#   pars = c('theta_y', 'theta_m', 'b0', 'b_size', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # mutiple years, different betas
# dat_stan$clim1_means <- rowMeans( clim_pred[,1:12] )
# dat_stan$clim2_means <- rowMeans( clim_pred[,13:24] )
# dat_stan$clim3_means <- rowMeans( clim_pred[,25:36] )
# 
# # Multiple years, weighted with simplex
# fit_year_simpl <- stan(
#   file = paste0("analysis/stan/grow/gaussian_mYears_simplex.stan"),
#   data = dat_stan,
#   pars = c('theta_k', 'b0', 'b_size', 'b_c', 's'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # Multiple years, weighted with simplex
# fit_year_simpl_re <- stan(
#   file = paste0("analysis/stan/grow/gaussian_mYears_simplex_re.stan"),
#   data = dat_stan,
#   pars = c('theta_k', 'b0', 's_yr', 'b_yr', 'b_size', 's', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # update data list
# dat_stan$K            <- 2
# dat_stan$M            <- 12
# 
# # power exponential moving window
# fit_24_nest <- stan(
#   file = paste0("analysis/stan/bernoulli_dirichlet_nest_24.stan"),
#   data = dat_stan,
#   pars = c('theta_y', 'theta_m', 'b0', 'b_c'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # average climate
# dat_stan$clim_means <- rowMeans( clim_pred )
# fit_avg <- stan(
#   file = paste0("analysis/bernoulli_avg_site.stan"),
#   data = dat_stan,
#   pars = c('beta', 'beta_site'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # Multiple years, weighted with simplex
# fit_mYear_simpl <- stan(
#   file = paste0("analysis/bernoulli_mYears_simplex_site.stan"),
#   data = dat_stan,
#   pars = c('theta', 'beta', 'beta_site'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# # average climate of first year
# dat_stan$clim_means <- rowMeans( clim_pred[,1:12] )
# fit_avg_yr1 <- stan(
#   file = paste0("analysis/stan/bernoulli_avg_site.stan"),
#   data = dat_stan,
#   pars = c('beta', 'beta_site', 'log_lik'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # average climate of first year
# dat_stan$clim_means <- rowMeans( clim_pred[,13:24] )
# fit_avg_yr2 <- stan(
#   file = paste0("analysis/stan/bernoulli_avg_site.stan"),
#   data = dat_stan,
#   pars = c('beta', 'beta_site', 'log_lik'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # average climate of second year
# dat_stan$clim_means <- rowMeans( clim_pred[,25:36] )
# fit_avg_yr3 <- stan(
#   file = paste0("analysis/stan/bernoulli_avg_site.stan"),
#   data = dat_stan,
#   pars = c('beta', 'beta_site', 'log_lik'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# 
# 
# # Multiple years, weighted with simplex
# fit_mYear_simpl <- stan(
#   file = paste0("analysis/bernoulli_mYears_simplex_site.stan"),
#   data = dat_stan,
#   pars = c('theta', 'beta', 'beta_site'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # power exponential moving window
# fit_gev <- stan(
#   file = paste0("analysis/bernoulli_gev_site.stan"),
#   data = dat_stan,
#   pars = c('loc', 'scale', 'shape', 'beta_site', 'beta'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# 
# dat_stan$clim <- t(clim_pred)
# # power exponential moving window
# fit_24 <- stan(
#   file = paste0("analysis/bernoulli_dirichlet_site.stan"),
#   data = dat_stan,
#   pars = c('theta', 'beta', 'beta_site'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # power exponential moving window
# dat_stan$clim <- clim_pred
# fit_expp <- stan(
#   file = paste0("analysis/bernoulli_expp.stan"),
#   data = dat_stan,
#   pars = c('sens_mu', 'sens_sd', 'beta_site', 'beta'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # gaussian moving window
# fit_gaus <- stan(
#   file = paste0("analysis/bernoulli_gaus.stan"),
#   data = dat_stan,
#   pars = c('sens_mu', 'sens_sd', 'alpha', 'beta'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# # # gaussian moving window
# # fit_null <- stan(
# #   file = paste0("analysis/bernoulli_null.stan"),
# #   data = dat_stan,
# #   pars = c('alpha'),
# #   warmup = sim_pars$warmup,
# #   iter = sim_pars$iter,
# #   thin = sim_pars$thin,
# #   chains = sim_pars$chains#,
# #   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# # )
# 
# 
# # plot it out -------------------------------------------------------------------------------------
# theta       <- extract(fit_24)[['theta']] %>% 
#                   as.data.frame %>% 
#                   stack %>% 
#                   mutate( ind = as.character(ind) ) %>%
#                   mutate( ind = gsub("V","", ind) ) %>%
#                   mutate( ind = as.numeric(ind) )
# theta_mean  <- summary(fit_24)$summary[,'mean'][paste0('theta[',1:24,']')]
# beta0_mean  <- summary(fit_24)$summary[,'mean'][paste0('beta_site[',1:7,']')] %>% mean
# 
# 
# xs          <- data.frame( year = c(2008:2016), 
#                            xs = rowSums( sweep(clim_pred, 2, theta_mean, '*') ) )
# 
# ys          <- seedl %>% 
#                   group_by(year, location) %>% 
#                   summarise( surv_sum = sum(surv_t1),
#                              rep      = n() ) %>%
#                   mutate( prop_surv = surv_sum / rep )
# xs      <- data.frame( year = c(2008:2016), x1 = rowMeans(clim_pred[,1:12]), 
#                                             x2 = rowMeans(clim_pred[,13:24]) )
# plot_d  <- full_join(ys, xs)
# x_seq   <- seq(min(plot_d$xs), max(plot_d$xs), length.out=100 )
# pred    <- boot::inv.logit( beta0_mean + x_seq * beta0_mean )
# 
# 
# 
# tiff("results/Seedling_survival_prec_dirichlet.tiff", unit="in", width=4.5, height=6.3, 
#      res=600,compression="lzw")
# 
# par(mfrow = c(2,1), mar = c(3,3,0.1,0.1), mgp = c(1.6,0.7,0) )
# boxplot(values ~ ind, data = theta, outline = F,
#         ylab = "Month weights", xlab = "Month before demographic observation",
#         cex.lab = 1)
# abline(v = 12.5, lty = 2)
# plot(prop_surv ~ xs, data = plot_d, col = as.factor(ys$location), pch = 16,
#      ylab = "Average seedling survival (year t)",  
#      xlab = "Average climate predictor", cex.lab = 1)
# lines(x_seq, pred, lwd = 2)
# 
# dev.off()
# 
# 
# 
# 
# x_seq1 <- seq(min(plot_d$x1), max(plot_d$x1), length.out=100)
# x_seq2 <- seq(min(plot_d$x2), max(plot_d$x2), length.out=100)
# beta0  <- summary(fit_mYear_beta)$summary[,'mean'][paste0('beta_site[',1:7,']')] %>% mean
# beta1  <- summary(fit_mYear_beta)$summary[,'mean']['beta[1]']
# beta2  <- summary(fit_mYear_beta)$summary[,'mean']['beta[2]']
# pred1  <- boot::inv.logit(beta0 + x_seq1 * beta1)
# pred2  <- boot::inv.logit(beta0 + x_seq2 * beta2)
# 
# # get posterior
# fit_extract <- extract(fit_mYear_beta)
# beta1_p     <- fit_extract$beta[,1]
# beta2_p     <- fit_extract$beta[,2]
# beta0_p     <- rowMeans(fit_extract$beta_site)
# pred1_p     <- lapply(1:6000, function(ii) boot::inv.logit(beta0_p[ii] + x_seq1 * beta1_p[ii])) %>%
#                   rbind_l
# pred2_p     <- lapply(1:6000, function(ii) boot::inv.logit(beta0_p[ii] + x_seq2 * beta2_p[ii])) %>%
#                   rbind_l
# 
# plot_post   <- function(ii, x_vals, pred_vec){
#   lines(x_vals, pred_vec[ii,], col="grey")
# }
#   
# 
# tiff("results/Seedling_survival_precipitation.tiff", unit="in", width=4.5, height=6.3, 
#      res=600,compression="lzw")
# 
# par(mfrow = c(2,1), mar = c(3,3,0.1,5), mgp = c(1.7,0.7,0) )
# plot(prop_surv ~ x1, data=plot_d, ylim = c(0,1), pch = 16, col = as.factor(ys$location),
#      ylab = "Average seedling survival (year t)", xlab = "Precipitation anomaly (year t)", type = "n")
# lapply(1:6000, plot_post, x_seq1, pred1_p)
# points(prop_surv ~ x1, data=plot_d, ylim = c(0,1), pch = 16, col = as.factor(ys$location),
#      ylab = "Average seedling survival (year t)", xlab = "Precipitation anomaly (year t)")
# lines( x_seq1, pred1, lwd=2 )
# 
# 
# legend(0.75,1, unique(ys$location), pch =16, col = unique(as.factor(ys$location)), bty = 'n',
#        xpd=T)
# plot(prop_surv ~ x2, data=plot_d, ylim = c(0,1), pch = 16, col = as.factor(ys$location),
#      ylab = "Average seedling survival (year t)", xlab = "Precipitation anomaly (year t-1)", type = "n")
# lapply(1:6000, plot_post, x_seq2, pred2_p)
# points(prop_surv ~ x2, data=plot_d, ylim = c(0,1), pch = 16, col = as.factor(ys$location),
#      ylab = "Average seedling survival (year t)", xlab = "Precipitation anomaly (year t-1)")
# lines( x_seq2, pred2, lwd=2 )
# legend(-0.2,0.2, c("posterior", "mean"), lwd = c(1,2),  col = c("grey","black"), bty = 'n',
#        xpd=T)
# 
# dev.off()
# 
# 
# 
# # parameter values and diagnostics ----------------------------------------------------------------
# 
# # list of model fits
# mod_fit   <- list(fit_1 = fit_site,
#                   fit_2 = fit_mYear_beta, 
#                   fit_2 = fit_avg_yr1,
#                   fit_4 = fit_avg_yr2)
# 
# # parameter values
# pars      <- c('beta', 'beta_site', 'log_lik')
# 
# # get central tendencies
# pars_diag_extract <- function(x){
#   
#   # central tendencies
#   tmp         <- rstan::extract(x)
#   par_means   <- sapply(tmp, function(x) mean(x)) %>%
#     setNames( paste0(names(tmp),"_mean") )
#   par_medians <- sapply(tmp, function(x) median(x)) %>%
#     setNames( paste0(names(tmp),"_median") )
#   central_tend<- c(par_means, par_medians)
#   
#   # diagnostics
#   diverg      <- do.call(rbind, args = get_sampler_params(x, inc_warmup = F))[,5]
#   n_diverg    <- length(which(diverg == 1))
#   df_summ     <- as.data.frame(summary(x)$summary)
#   rhat_high   <- length(which(df_summ$Rhat > 1.1))
#   n_eff       <- df_summ$n_eff / length(diverg)
#   n_eff_low   <- length(which(n_eff < 0.1))
#   mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
#   diagnostics <- c(n_diverg = n_diverg, rhat_high = rhat_high,
#                    n_eff_low = n_eff_low, mcse_high = mcse_high)
#   out         <- c( central_tend, diagnostics ) %>% t %>% as.data.frame
#   
#   rm(tmp) ; return(out)
#   
# }
# 
# # store posteriors
# posterior_extract <- function(model_fit, model_name){
#   
#   # central tendencies
#   tmp         <- rstan::extract(model_fit)
#   post_df     <- do.call(cbind, tmp) %>% as.data.frame
#   ll_id       <- grep("V", colnames(post_df) )
#   new_names   <- paste0("log_lik_", 1:length(ll_id) )
#   names(post_df)[ll_id] <- new_names # no way to do this in dplyr
#   post_df     <- tibble::add_column(post_df,
#                                     model = model_name, .before=1)
#   
#   rm(tmp) ; return(post_df)
#   
# }
# 
# # calculate central tendencies
# pars_diag_l   <- lapply(mod_fit, pars_diag_extract)
# mod_pars_diag <- Reduce(function(...) bind_rows(...), pars_diag_l) %>%
#                     tibble::add_column(model = names(mod_fit), .before = 1)
# 
# # store posteriors
# posts_l       <- Map(posterior_extract, mod_fit, names(mod_fit) )
# posteriors    <- Reduce(function(...) bind_rows(...), posts_l)
# 
# 
# # WAIC model comparison --------------------------------------------------------------------
# 
# # wAIC model selection using loo approximation (from library 'loo')
# log_liks  <- lapply(mod_fit, extract_log_lik)
# 
# # leave-one-out estimates
# loo_l      <- lapply(log_liks, loo) %>%
#                 setNames( c('loo_fit1',  'loo_fit2', 'loo_fit3',  'loo_fit4') )
# loo_df     <- loo::compare(loo_l$loo_fit1, loo_l$loo_fit2, loo_l$loo_fit3, loo_l$loo_fit4) %>%
#                 as.data.frame %>%
#                 tibble::add_column(model = gsub("loo_","",names(loo_l) ), .before = 1)
# 
# # WAIC estimates
# waic_l    <- lapply(log_liks, waic) %>%
#                 setNames( c('waic_fit1',  'waic_fit2', 'waic_fit3',  'waic_fit4') )
# waic_df   <- loo::compare(waic_l$waic_fit1, waic_l$waic_fit2, waic_l$waic_fit3, waic_l$waic_fit4) %>%
#                 as.data.frame %>%
#                 tibble::add_column(model = gsub("waic_","",names(waic_l) ), .before = 1)
# 
