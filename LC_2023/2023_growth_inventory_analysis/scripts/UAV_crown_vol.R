# Process non-thermal UAV data 

library(ggplot2)
library(viridis)
library(dplyr)

#import data from UAV flight in August 2022. Prepared by Cory Garms

UAV_Data <- read.csv(file = "LC_2022/2022_Marchel_drone/2022_UAV_tables/LC_full_complete.csv")

UAV_indices <- UAV_Data[c(2:44)]




#compare UAV crown index derived from shadow sizes to tree volume.
#interested in V497 for this comparison, volume near to when drone imagery taken.

growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")
UAV_join <- inner_join(UAV_Data, growth_2022, by = "ID")

growth_2022 <- subset(growth_dat, select = c("ID","event_short","V497","H497","D497","block"))
growth_2022$HD <- growth_2022$H497/growth_2022$D497

UAV_join <- inner_join(UAV_Data, growth_2022, by = "ID")

#scatter plot

#log transformed
ggplot(UAV_join, aes(x=log(V497),y=log(shadow_area), color = HD))+
  geom_point()+
  scale_color_viridis()+
  geom_abline(intercept=-4.2078,slope = 0.4972)
  

#would expect that trees with a lower height to diameter ratio may be branchier and tend to be above the fitted line.
?geom_abline

branchiness_lm <- lm(log(shadow_area)~log(V497), data =UAV_join )
summary(branchiness_lm)
plot(branchiness_lm, which = 1)
#outliers to look at to see if they are branchier than expected

#residuals from this relationship could feasibly be related to deviations of crown area from stem volume due to branchiness

ggplot(UAV_join, aes(y=HD, x=branchiness_lm$residuals))+
         geom_point()

#there may be a slight relationship here but HD ratio does not seem to explain branchiness residuals

#how to branchiness residuals vary spatially?

library(latticeExtra)

levelplot(branchiness_lm$residuals ~-column*row, UAV_join,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
#this may be helpful. It does seem to match with where I see branchier trees.

levelplot(shadow_area ~-column*row, UAV_join,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))


library(ggmap)
library(terra)
library(sf)
nc <- read_sf(system.file("gpkg/nc.gpkg", package = "sf"))
plot(nc["BIR74"], reset = FALSE, key.pos = 4)
plot(st_buffer(nc[1,1], units::set_units(10, km)), col = 'NA', 
     border = 'red', lwd = 2, add = TRUE)

install.packages("tidyverse")
library(tidyverse)

nc.32119 <- st_transform(nc, 32119) 
year_labels <- 
  c("SID74" = "1974 - 1978", "SID79" = "1979 - 1984")
nc.32119 |> select(SID74, SID79) |> 
  pivot_longer(starts_with("SID")) -> nc_longer

ggplot() + geom_sf(data = nc_longer, aes(fill = value), linewidth = 0.4) + 
  facet_wrap(~ name, ncol = 1, 
             labeller = labeller(name = year_labels)) +
  scale_y_continuous(breaks = 34:36) +
  scale_fill_gradientn(colours = sf.colors(20)) +
  theme(panel.grid.major = element_line(colour = "white"))


str(UAV_join)
ggplot() + geom_sf(data = UAV_join)
SpatialPoints()