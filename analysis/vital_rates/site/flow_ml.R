# simplified script to test for SITE BY CLIMATE interaction
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
flow <- subset(lupine_df, !is.na(flow_t1) ) %>% 
          subset( area_t1 != 0) %>% 
          mutate( log_area_t1  = log(area_t1),
                  log_area_t12 = log(area_t1)^2,
                  year         = year + 1 ) #%>% 
          #mutate( log_area_t1  = scale(log_area_t1) %>% as.vector,
          #        log_area_t12 = scale(log_area_t12) %>% as.vector )

        
# climate format ----------------------------------------------------------------
years     <- unique(flow$year)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1','_t0_tm1','_t0_tm2') )
  
  mutate(x, 
         avgt0     = x %>% select(V1:V12) %>% rowSums,
         avgtm1    = x %>% select(V13:V24) %>% rowSums,
         avgt0_tm1 = x %>% select(V1:V24) %>% rowSums,
         avgt0_tm2 = x %>% select(V1:V36) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1, avgt0_tm1, avgt0_tm2) %>% 
    setNames( c('year',clim_names) )
  
}

# format climate - need to select climate predictor first 
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs) %>% 
              year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp')
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_anom('oni')

spei_mat <- subset(clim, clim_var == 'spei' ) %>%
              spei_clim_form(years, m_back, m_obs) %>% 
              year_anom('spei')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat,spei_mat) )

# demography plus clim
flow_clim  <- left_join(flow, clim_mat) %>%
                subset( !is.na(location) )



# fit climate models -------------------------------------------------
climate_mods <- list(

  # null 
  flow_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location),
  
  # precipitation
  # flow_t1 ~ log_area_t1 + log_area_t12 + ppt_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + ppt_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + ppt_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + ppt_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # temperature
  flow_t1 ~ log_area_t1 + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location),
  flow_t1 ~ log_area_t1 + tmp_tm1 + (log_area_t1 | year) + (log_area_t1 | location),
  flow_t1 ~ log_area_t1 + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location) +
            (tmp_t0 | location),
  flow_t1 ~ log_area_t1 + tmp_tm1 + (log_area_t1 | year) + (log_area_t1 | location) +
            (tmp_tm1 | location)
  # flow_t1 ~ log_area_t1 + log_area_t12 + tmp_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + tmp_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # # enso
  # flow_t1 ~ log_area_t1 + log_area_t12 + oni_t0 + (log_area_t1 | year) + (log_area_t1 | location)+ (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + oni_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + oni_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + oni_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # # spei
  # flow_t1 ~ log_area_t1 + log_area_t12 + spei_t0 + (log_area_t1 | year) + (log_area_t1 | location)+ (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + spei_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + spei_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  # flow_t1 ~ log_area_t1 + log_area_t12 + spei_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid)

)


# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data = flow_clim, family='binomial') ) %>% 
              setNames( c( 'null',
                           'tmp_t0',  'tmp_tm1',  'tmp_t0_site',  'tmp_tm1_site') )

AICtab(mod_clim, weights=T)

# out data frame
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/flow/flow_mod_sel.csv', row.names = F)


# tmp_t0 
best_mod <- glmer(flow_t0 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location), 
                  data = flow_clim, family='binomial') 
#best_mod <- mod_clim[[out$model[1]]]


# random effects
re_df   <- ranef(best_mod)$year %>% 
              tibble::add_column(.before=1, coef = row.names(.) ) %>% 
              bind_rows( ranef(best_mod)$location %>% 
                         tibble::add_column(.before=1, coef = row.names(.) )
              ) %>% 
              # bind_rows( ranef(best_mod)$newid %>% 
              #            tibble::add_column(.before=1, coef = row.names(.) )
              # ) %>% 
              mutate( type_coef = 'ranef' )

# out (binding fixed effects)
out_df  <- fixef(best_mod) %>% 
              t %>% t %>% 
              as.data.frame %>% 
              tibble::add_column(.before=1, ranef = row.names(.) ) %>% 
              mutate( type_coef = 'fixef' ) %>% 
              bind_rows( re_df )
     
write.csv(out_df, 
          'results/ml_mod_sel/flow/flow_best_mod.csv',
          row.names=F)


# plot -----------------------------------------------------------
coefs       <- best_mod %>% fixef

# quantiles of climate predictor
clim_quant  <- flow_clim$tmp_t0 %>% unique %>% quantile

x_seq <- seq( min(flow_clim$log_area_t1, na.rm=T),
              max(flow_clim$log_area_t1, na.rm=T),
              length.out = 100 )

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['25%']

y_high <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['75%']


tiff('results/ml_mod_sel/flow.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(jitter(flow_t1) ~ log_area_t0, data = flow_clim,
     ylab = 'survival probability')
lines(x_seq, boot::inv.logit(y_low), lwd=2, col = 'blue')
lines(x_seq, boot::inv.logit(y_high), lwd=2, col = 'red')

legend(-0.5,0.8, 
       c('low tmp t0','high tmp t0'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()


# exploratory graphs --------------------------------------------------------------------

# sample sizes
rep_df  <- flow_clim %>% 
              group_by( year,location ) %>% 
              summarise( succ  = sum(flow_t1, na.rm=T),
                         fail  = sum(flow_t1 == 0, na.rm=T) ) %>% 
              ungroup %>% 
              mutate( tot = succ + fail,
                      n_i = round(succ/(succ+fail),1),
                      y = 0.8) %>% 
              as.data.frame %>% 
              select(year,location, n_i, tot, y ) #%>% 


# Make the multi-panel plot 
ggplot(data  = flow_clim, 
       aes(x = log_area_t1, 
           y = flow_t1) ) +
  geom_point(alpha = 0.5,
             pch = 1) +
  # Fit linear Regression
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"),
              size = 0.5) +
  # split in panels
  facet_grid(location ~ year) +
  theme_bw() +
  ggtitle("Flowering (flow_t1 <> NA, log_area_t1 <> 0)") +
  geom_text(
    data    = rep_df,
    mapping = aes(x = -Inf, y = -Inf, label = n_i),
    hjust   = -0.1,
    vjust   = -1
  ) +
   geom_text(
    data    = rep_df,
    mapping = aes(x = -Inf, y = y, label = tot),
    hjust   = -0.1,
    vjust   = -1
  ) +
  ggsave(filename = "results/ml_mod_sel/flow_checks.png",
         dpi = 300, width = 50, height = 30, units = "cm")



# year pooled data checks -----------------------------------------------------
pooled_df  <- flow_clim %>% 
                  group_by( year ) %>% 
                  summarise( succ  = sum(flow_t1, na.rm=T),
                             fail  = sum(flow_t1 == 0, na.rm=T) ) %>% 
                  ungroup %>% 
                  mutate( tot  = succ + fail,
                          prop = round(succ/(succ+fail),3) ) %>% 
                  as.data.frame %>% 
                  select(year, prop )

# pool climate and survival info together
pool_df <- left_join(pooled_df,clim_mat) 

    
# plot by variable
plot_by_var <- function(var){
  
  # nonstandard evaluation
  var2 <- rlang::enquo(var)
  
  ggplot(data  = pool_df, 
         aes(x = !! var2, 
             y = prop,
             size = 2) ) +
    ylab('Flowering probability')+
    theme(axis.title=element_text(size=20),
          legend.position="none") +
    geom_point(alpha = 1,
               pch   = 16)
}
   
# plot it all out 
library(gridExtra)
p <- grid.arrange(plot_by_var(ppt_t0), plot_by_var(ppt_tm1), 
                  plot_by_var(ppt_t0_tm1), plot_by_var(ppt_t0_tm2), 
                  plot_by_var(tmp_t0), plot_by_var(tmp_tm1), 
                  plot_by_var(tmp_t0_tm1), plot_by_var(tmp_t0_tm2), 
                  plot_by_var(oni_t0), plot_by_var(oni_tm1), 
                  plot_by_var(oni_t0_tm1), plot_by_var(oni_t0_tm2),
                  nrow = 3, ncol=4)

ggsave(filename = "results/ml_mod_sel/flow_tab.png",
       p, 
       dpi = 300, width = 50, height = 30, units = "cm")
