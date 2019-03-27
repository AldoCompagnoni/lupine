rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(testthat)

clim_var <- 'oni'

result_f <- list.files('results') %>% 
              grep('.csv',.,value=T) %>% 
              grep(clim_var,.,value=T,perl=T) %>% 
              grep('^((?!_12m_).)*$',.,value=T,perl=T)  

result_f_post <- grep('_post_',result_f,value=T,perl=T)
result_f_summ <- grep('_summ_',result_f,value=T,perl=T)

# # no RE
# avg_f <- result_f %>% 
#             grep('prec',.,value=T) %>% 
#             grep('_a_',.,value=T) %>% 
#             grep('^((?!_re_).)*$',.,value=T,perl=T) %>% 
#             grep('^((?!_post_).)*$',.,value=T,perl=T)  
# 
# y3_f  <- result_f %>% 
#             grep('prec',.,value=T) %>% 
#             grep('_3y_',.,value=T) %>% 
#             grep('^((?!_re_).)*$',.,value=T,perl=T) %>% 
#             grep('^((?!_post_).)*$',.,value=T,perl=T)  

# Random effect 
avg_re_f <- result_f %>% 
              grep(clim_var,.,value=T) %>% 
              grep('_a_',.,value=T) %>% 
              grep('_re_',.,value=T,perl=T) %>% 
              grep('^((?!_post_).)*$',.,value=T,perl=T)  
              
y3_re_f  <- result_f %>% 
              grep(clim_var,.,value=T) %>% 
              grep('_3y_',.,value=T) %>% 
              grep('_re_',.,value=T,perl=T) %>% 
              grep('^((?!_post_).)*$',.,value=T,perl=T)  


# Compare coeffs -------------------------------------------------

# check convergence by reading all summaries
all_summ <- grep('^((?!_post_).)*$',result_f,value=T,perl=T) %>% 
              paste0('results/',.) %>% 
              lapply(read.csv) %>% 
              lapply(function(x) select(x,Rhat)) %>% 
              lapply(function(x) any( x > 1.05 )) %>% 
              unlist

conv_prob<- grep('^((?!_post_).)*$',result_f,value=T,perl=T) %>% 
              .[which(all_summ)]


# get to actual coefficients -----------------------------------------------

# avg_l     <- paste0('results/',avg_f) %>% 
#               lapply( read.csv ) %>% 
#               lapply(function(x) subset(x, coeff == 'b_c') %>% 
#                                  select(coeff,mean,sd) ) %>% 
#               setNames( gsub('_a_summ_prec.csv','',avg_f) )
#         
avg_re_l  <- paste0('results/',avg_re_f) %>%
              lapply( read.csv ) %>%
              lapply(function(x) subset(x, X == 'b_c') %>%
                                 select(X,mean,sd) ) %>%
              setNames( gsub('_a_summ_re_prec.csv','',avg_re_f) )
              
# three years weighted
# y3_l     <- paste0('results/',y3_f) %>% 
#               lapply( read.csv ) %>% 
#               lapply(function(x) subset(x, grepl('theta',coeff) | coeff == 'b_c') %>% 
#                                  select(coeff,mean,sd) ) %>% 
#               setNames( gsub('_sam_3y_summ_prec.csv','',y3_f) )
y3_re_l  <- paste0('results/',y3_re_f) %>% 
              lapply( read.csv ) %>% 
              lapply(function(x) subset(x, grepl('theta',X) | X == 'b_c') %>% 
                                 select(X,mean,sd) ) %>% 
              setNames( gsub('_sam_3y_summ_re_airt.csv','',y3_re_f) )



# plots B_C posterior ----------------------------------------------------------

# 3-year
posts_vec <- result_f_post %>% 
                grep('_3y_',.,value=T) %>% 
                grep('_re_',.,value=T)  
                
b_c_post_3y <- lapply(posts_vec, function(x) read.csv(paste0('results/',x))) %>% 
                  lapply(function(x) dplyr::select(x,b_c)) %>% 
                  do.call(cbind, .) %>% 
                  # alphabetical order
                  setNames( gsub(paste0('_sam_3y_post_re_',clim_var,'.csv'),'',
                                 posts_vec) ) %>% 
                  stack

# averages
posts_vec <- result_f_post %>% 
                grep('_a_',.,value=T) %>% 
                grep('_re_',.,value=T)  
                
b_c_4       <- lapply(posts_vec[-4], function(x) read.csv(paste0('results/',x))) %>% 
                lapply(function(x) dplyr::select(x,b_c)) %>% 
                do.call(cbind, .) %>% 
                # alphabetical order
                setNames( gsub(paste0('_a_post_re_',clim_var,'.csv'),'',posts_vec[-4]) ) %>% 
                stack %>% 
                mutate( ind = as.character(ind) )
b_c_post_a  <-  lapply(posts_vec[4], function(x) read.csv(paste0('results/',x))) %>% 
                lapply(function(x) dplyr::select(x,b_c)) %>%
                .[[1]] %>% 
                mutate( ind = 'sl' ) %>% 
                setNames( c('values','ind') ) %>% 
                rbind( b_c_4 )
                

tiff(paste0('results/figs/',clim_var,'_climate_coeffs.tiff'),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

# plot out
par(mfrow=c(2,1), mar=c(2.8,4,0.1,0.1))
boxplot(values ~ ind, data = b_c_post_3y,
        ylab = paste0(clim_var,' effect, 3-years weighted') )
abline(h=0, lty=2)

boxplot(values ~ ind, data = b_c_post_a,
        ylab = paste0(clim_var,' effect, 3-year average') )
abline(h=0, lty=2)

dev.off()


# 3-year
posts_vec <- result_f_post %>% 
                grep('_3y_',.,value=T) %>% 
                grep('_re_',.,value=T)  
                
tk_post_3y <- lapply(posts_vec, function(x) read.csv(paste0('results/',x))) %>% 
                  lapply(function(x) dplyr::select(x,theta_k.1:theta_k.3))


tk_df <- Map(function(x,y) mutate(x, vr = y), tk_post_3y,
             gsub('_sam_3y_post_re_oni.csv','',posts_vec) ) %>% 
              Reduce(function(...) rbind(...), .) %>% 
              gather(theta,value, -vr)

tiff(paste0("results/figs/theta_coeffs.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2))
boxplot(value ~ theta + vr, data= subset(tk_df, vr =='fert'))
boxplot(value ~ theta + vr, data= subset(tk_df, vr =='flow'))
boxplot(value ~ theta + vr, data= subset(tk_df, vr =='grow'))
boxplot(value ~ theta + vr, data= subset(tk_df, vr =='surv'))

dev.off()


b_c_post_3y <- lapply(b_c_post_3y, function(x))
                  do.call(cbind, .) %>% 
                  # alphabetical order
                  setNames( gsub('_sam_3y_post_re_prec.csv','',posts_vec) ) %>% 
                  stack
