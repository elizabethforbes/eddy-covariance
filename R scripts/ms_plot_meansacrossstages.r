# plotting the main results: mean fluxes by crop stage and management, only
# the empirical data (no gap-filled days):
library(tidyverse)
library(patchwork)

source(here::here("R scripts", "plot_databymanagement.r"))

# plot NEE
p2 <- plot_management_comparison(data1, flux_var = "NEE")

# plot GPP
g2 <- plot_management_comparison(data1, flux_var = "GPP")

# plot Reco
r2 <- plot_management_comparison(data1, flux_var = "Reco")
r2 + labs( 
  y = (bquote(respiration ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")),
  title = "Respiration by Management and Crop Stage")

# plot all three: remove titles, merge legends
p2_a <- p2 + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                   legend.position = "blank")
g2_a <- g2 + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                   legend.position = "blank")
r2_a <- r2 + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                   legend.position = "bottom")

p2_a/g2_a/r2_a

# save with width = 500, height = 900