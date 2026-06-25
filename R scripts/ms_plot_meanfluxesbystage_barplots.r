source(here::here("R scripts", "plot_databymanagement.r"))

library(tidyverse)
library(patchwork)

n <- plot_management_comparison(data1, flux_var = "NEE")
g <- plot_management_comparison(data1, flux_var = "GPP")
r <- plot_management_comparison(data1, flux_var = "Reco")+labs( 
  y = (bquote(respiration ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")),
  title = "Respiration by Management and Crop Stage")
n
g
r

# plot together:
n2 <- n + labs(title = element_blank()) +
  theme(legend.position = "none")
g2 <- g +
  labs(title = element_blank(), subtitle = element_blank()) +
  theme(legend.position = "none")
r2 <- r +
  labs(title = element_blank(), subtitle = element_blank())+
  theme(legend.position = "bottom")

n2/g2/r2
