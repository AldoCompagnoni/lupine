# /////////////////////////////////////////////////////////////////////////
#
# Explorative graphs for some climate variables: ppt, tmp, oni and tmp across
# time.
#
# Compile script as report CTRL+SHIFT+K and move results in:
# ~\Dropbox\lupine\results\climate
#
# Valentin Stefan
#
# /////////////////////////////////////////////////////////////////////////

library(GGally)
library(plotly)
library(data.table)

rm(list = ls(all.names = TRUE))

# The hack with the setting of the working directory below is important only if
# compiling this script as HTML report when running from within RStudio (via
# lupine.Rproj).
setwd(gsub(pattern = '/analysis/get_climate', replacement = '', x = getwd()))

dt <- read.csv("climate_demographic_year_lupine.csv")


# Matrix plots ------------------------------------------------------------

# Pairs

pairs <- ggpairs(dt[, c("ppt", "tmp", "oni")],
                 aes(text = dt$year),
                 lower = list(continuous = "smooth"))
ggplotly(pairs)
ggsave(filename = "results/climate/pairs_ppt_tmp_oni.png", 
       plot = pairs, 
       width = 20, height = 15, units = "cm", dpi = 300)


# Individual plots

p1 <- ggplot(dt, aes(x = ppt, y = tmp, text = year)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw()
ggplotly(p1)
ggsave(filename = "results/climate/pairs_ppt_vs_tmp.png", 
       plot = p1, 
       width = 10, height = 8, units = "cm", dpi = 300)

p2 <- ggplot(dt, aes(x = ppt, y = oni, text = year)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw()
ggplotly(p2)
ggsave(filename = "results/climate/pairs_ppt_vs_oni.png", 
       plot = p2, 
       width = 10, height = 8, units = "cm", dpi = 300)

p3 <- ggplot(dt, aes(x = tmp, y = oni, text = year)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw()
ggplotly(p3)
ggsave(filename = "results/climate/pairs_tmp_vs_oni.png", 
       plot = p3, 
       width = 10, height = 8, units = "cm", dpi = 300)


# ppt, tmp, oni & anomalies across time  ----------------------------------

dt$yr_lbl <- paste(dt$year - 1, dt$year, sep = "-")
dt$year <- NULL
dtm <- data.table::melt(dt, id.vars = "yr_lbl")

dtm$group <- factor(dtm$variable, 
                    levels = c("ppt", "ppt_anom", "tmp", "tmp_anom", "oni", "oni_anom"),
                    labels = rep(c("real", "anomaly"), times = 3),
                    ordered = TRUE)
dtm$var <- factor(dtm$variable, 
                  levels = c("ppt", "ppt_anom", "tmp", "tmp_anom", "oni", "oni_anom"),
                  labels = rep(c("ppt", "tmp", "oni"), each = 2),
                  ordered = TRUE)

time_series <- ggplot(dtm, aes(x = yr_lbl, y = value)) +
  geom_point() +
  geom_line(group = 1) +
  facet_wrap(vars(var, group), ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplotly(time_series)

ggsave(filename = "results/climate/time_series_ppt_tmp_oni_and_anomalies.png", 
       plot = time_series, 
       width = 21, height = 15, units = "cm", dpi = 300)


# Copy the html report into '/results/climate/'
file.copy(from = paste0(getwd(), '/analysis/get_climate/explorative_clim.html'),
          to = paste0(getwd(), '/results/climate/explorative_clim.html'),
          overwrite = TRUE)
