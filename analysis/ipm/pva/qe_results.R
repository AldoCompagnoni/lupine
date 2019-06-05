# graphs comparison between quasi-extinction sims /w or w/o demographic stochasticity
library(dplyr)
library(tidyr)
library(ggplot2)
options(stringsAsFactors = F)

# file paths
f_pth  <- paste0( 'results/ipm/qe_sims/eve/',
                  list.files('results/ipm/qe_sims/eve') )
                  
# put together
all_ds_sim <- lapply( f_pth, read.csv) %>% 
              bind_rows %>% 
              arrange(sim_i, clim_x, loc_v, time)

# quasi-extinction with demographic stochasticity 
qe_ds_df <- all_ds_sim %>% 
                group_by(sim_i,clim_x,loc_v) %>% 
                summarise( qe = any( n < 5) ) %>% 
                ungroup %>% 
                mutate( qe = as.numeric(qe) ) %>% 
                group_by( clim_x,loc_v ) %>% 
                summarise( qe_p = sum(qe) / 1000 ) %>% 
                ungroup %>% 
                rename( Location = loc_v,
                        qe_ds_p  = qe_p )

# load QE results from no-demographic stochasticity simulations
load( 'results/ipm/qe_sims/qe_5.RData' )

# calculate quasi extinction probability 
qe_df <- mutate(sim_df, qe = unlist(qe_l) ) %>% 
            group_by( clim_x, loc_v ) %>% 
            summarise( qe_p = sum(qe) ) %>% 
            ungroup %>% 
            mutate( qe_p = qe_p / 1000 ) %>% 
            rename( Location    = loc_v )

# PVA with demographic stocasticity 
pva_ds_df  <- full_join( pva_df, qe_ds_df ) %>% 
                gather( sim_type, qe, qe_p:qe_ds_p ) %>% 
                mutate( sim_type = replace(sim_type, 
                                           sim_type == 'qe_ds_p',
                                           '0: demo stch.') ) %>% 
                mutate( sim_type = replace(sim_type, 
                                           sim_type == 'qe_p',
                                           '2: no demo stch.') )

pva_ds_df$sim_type %>% unique


# quasi-extinction probabilities --------------------------
ggplot(pva_df) +
  geom_line( aes(x     = clim_x,
                 y     = qe_p,
                 color = Location),
             size = 1,
             lty  = 5 ) + 
  viridis::scale_color_viridis( discrete=T ) +
  geom_line( data = qe_ds_df,
             aes(x     = clim_x,
                 y     = qe_ds_p,
                 color = Location),
             size = 1.5 ) +
  labs( x = 'Climate anomaly',
        y = 'Quasi-extinction probability') +
  ggsave( 'results/ipm/qe_sims/qe_5_demo_nodemo.tiff',
          width=6.3, height=6.3, compression='lzw' )
  

# single runs ----------------------------------------

# select run for DR (3) with certain climate
ids  <- sim_df$loc_v == 'DR (3)' & (sim_df$clim_x > -0.5323705 & sim_df$clim_x < 0) 

# get specific projections
dr_ss_l <- purrr::keep(proj_l, ids) 
dr_n_l  <- list()

for(ii in 1:length(dr_ss_l) ){
  
  dr_n_l[[ii]] <- 
    dr_ss_l[ii][[1]][3:102,] %>% colSums
  
}
  
# no demographic stochasticity   
dr_df <- dr_n_l %>% 
           bind_cols %>% 
           mutate( time = 1:nrow(.) ) %>% 
           gather( sim_i, n, V1:V1000) %>% 
           subset( time < 15 ) %>%
           mutate( sim_i = gsub('V','',sim_i) ) %>% 
           mutate( sim_i = as.integer(sim_i) ) %>% 
           mutate( sim   = 'no_dem_stoch' )
           

# with demographic stochasticity
dr_ds <- all_ds_sim %>% 
            subset(loc_v == 'DR (3)' & 
                     (clim_x > -0.5323705 & clim_x < 0) ) %>% 
            subset( time < 15 ) %>%
            select( sim_i, time, n ) %>% 
            mutate( sim   = 'no_dem_stoch' )


plot_df <- bind_rows(dr_df, dr_ds)

ggplot(plot_df) +
  geom_line( aes( x     = time,
                  y     = n,
                  group = sim_i,
                  color = sim,
                  ),
             alpha = 0.1 )



ggplot(dr_df) +
  geom_line( aes( x     = time,
                  y     = n,
                  group = sim_i),
             alpha = 0.1 ) +
  geom_line( data = dr_ds,
             aes( x     = time,
                  y     = n,
                  group = sim_i),
             alpha = 0.1,
             color = 'red' ) +
  labs( x = 'Time step',
        y = 'Population size' ) + 
  ggtitle( 'DR(3): Demographic stoch vs. no demo stoch') +
  ggsave( 'results/ipm/qe_sims/sim_DR(3).tiff',
          width=6.3, height=6.3, compression = 'lzw')









# plot 
ggplot(pva_ds_df) +
  geom_line( aes(x        = clim_x,
                 y        = qe,
                 color    = Location,
                 linetype = 'solid') ) #+
  scale_linetype_manual(values=c("longdash", "solid") ) 

# plot 
ggplot(pva_ds_df) +
  geom_line( aes(x     = clim_x,
                 y     = qe,
                 color = Location,
                 lty   = sim_type),
             size = 1.5) + 
  viridis::scale_color_viridis( discrete=T ) +
  scale_linetype_manual(values=c("solid", "longdash"))
