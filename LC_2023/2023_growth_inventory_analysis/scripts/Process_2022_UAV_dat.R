# Process UAV data from 2022

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)

#import data from UAV flight in August 2022. Prepared by Cory Garms

UAV_Data <- read.csv(file = "LC_2022/2022_Marchel_drone/2022_UAV_tables/LC_full_complete.csv")

UAV_indices <- UAV_Data[c(2:44)]
write.csv(UAV_indices,file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_indices.csv")
UAV_thermal <- UAV_Data[c(1:8,47:56)]
write.csv(UAV_thermal, "LC_2023/2023_growth_inventory_analysis/processed/UAV_thermal.csv")
UAV_allometry <- UAV_Data[c(1:8,57:59)]
write.csv(UAV_allometry, file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_allometry_2022.csv")

#Indices

#examine indices for outliers, correlations, imputtation

ggplot(data = UAV_indices) +
  geom_histogram(aes(x = TGI_mean))
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = (GRVI_mean)))
#skewed
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = (MGRVI_mean)))
#skewed
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = EXG_mean))
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = EXGR_mean))
#little skewed
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = NDVI_mean))
#little skewed
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = NDRE_mean))
ggplot(data = subset(UAV_indices, GNDVI_mean > 0.01)) +
  geom_histogram(aes(x = GNDVI_mean))

#possible outlier
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = SAVI_mean))
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = OSAVI_mean))
#skewed
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = MSAVI_mean))
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = GCI_mean))
ggplot(data = UAV_indices) +
  geom_histogram(aes(x = RECI_mean))


#examining relationships between indices
install.packages("cluster")
install.packages("factoextra")
library(cluster)
library(factoextra)
library(corrplot)

cor_set_indices <- subset(UAV_indices, select = c(ID,TGI_mean,GRVI_mean,MGRVI_mean,EXG_mean,EXGR_mean,NDVI_mean,NDRE_mean,GNDVI_mean,SAVI_mean,OSAVI_mean,MSAVI_mean,GCI_mean,RECI_mean))
colnames(cor_set_indices) <- c("ID","TGI","GRVI","MGRVI","EXG","EXGR","NDVI","NDRE","GNDVI","SAVI","OSAVI","MSAVI","GCI","RECI")

#simple correlation plot
corplt <- corrplot.mixed(cor(cor_set_indices[,-1]),upper = "shade",lower = "number")
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/figures/veg_indices_corr.png",plot=corplt,width = 8, height = 8, units = "in", dpi=300)


####heatmap

#scale
cor_set_indices_scaled <- as.data.frame(scale(cor_set_indices[,-1]))
cor_set_indices_scaled$ID <- cor_set_indices$ID
#reorder
cor_set_indices_scaled <- cor_set_indices_scaled[, c(14,1:13)]

#correlation matrix for heatmap
indices_cor <- cor((cor_set_indices_scaled[,-1]), use = "pairwise.complete.obs", method = "pearson")
#visualize

col<- colorRampPalette(c("red4", "white", "blue4"))(100)
pdf("LC_2023/2023_growth_inventory_analysis/figures/indices_heatmap.pdf")
indices_heatmap <-  heatmap(indices_cor, col=col, symm=TRUE)
dev.off()

### see if any of these indices explain variation in volume 
growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")
str(SPAD_dat)

#interest in volume at time of measurement and final volume
growth_UAV_comp <- subset(growth_dat, select = c("ID","V497","V801"))

indices_list <- list(growth_UAV_comp,cor_set_indices_scaled)
indices_growth <- indices_list %>% reduce(inner_join, by ="ID")


#simple corr plot
png(filename = "LC_2023/2023_growth_inventory_analysis/figures/cor_indices_growth.png",height = "8",width = "8",units = "in",res=300)
corrplot.mixed(cor(indices_growth[,-1]),upper = "shade",lower = "number")
dev.off()
#heat map 

growth_indices_cor <- cor((indices_growth[,-1]), use = "pairwise.complete.obs", method = "pearson")

#visualize

col<- colorRampPalette(c("red4", "white", "blue4"))(100)
pdf("LC_2023/2023_growth_inventory_analysis/figures/indices_growth_heatmap.pdf")
indices_heatmap_growth <-  heatmap(growth_indices_cor, col=col, symm=TRUE)
dev.off()

indices_heatmap_growth
growth_indices_cor[1:2,]

corrplot(growth_indices_cor[1:2,],col=col)


#####UAV allometry dat

#compare UAV crown index derived from shadow sizes to tree volume.
#interested in V497 for this comparison, volume near to when drone imagery taken.

growth_allometry <- subset(growth_dat, select = c("ID","event_short","V497","H497","D497","V801","H801","D801","block"))
growth_allometry$HD497 <- growth_allometry$H497/growth_allometry$D497
growth_allometry$HD801 <- growth_allometry$H801/growth_allometry$D801

UAV_allometry_join <- inner_join(UAV_allometry, growth_allometry, by = "ID")

#scatter plot
#compare crown area to height.
#trees that have larger crowns than predicted by height, may be expected to be branchier.
#would expect that trees with a lower height to diameter ratio may be branchier and tend to be above the fitted line.


#log transformed
branchiness_scatter <- ggplot(UAV_allometry_join, aes(x=log(shadow_area),y=log(H497), color = HD497))+
  geom_point()+
  scale_color_viridis()+
  geom_abline(intercept=8.34652,slope = 0.41282)

branchiness_scatter
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/figures/branchiness_scatter.png",branchiness_scatter,dpi = 300)
#larger trees that are below the line of fit have lower HD ratio -likely branchier
  
branchiness_lm_2022 <- lm(log(H497)~log(shadow_area), data =UAV_allometry_join )
summary(branchiness_lm_2022)
plot(branchiness_lm_2022, which = 1)

#outliers to look at to see if they are branchier than expected

#residuals from this relationship could feasibly be related to deviations of crown area from tree height due to branchiness

ggplot(UAV_allometry_join, aes(y=H801, x=branchiness_lm_2022$residuals))+
         geom_point()+
  xlab("branchiness resisdual")+
  ylab("final tree height (cm)")+
  theme_bw()


#how to branchiness residuals vary spatially?

library(latticeExtra)

levelplot(branchiness_lm_2022$residuals ~-column*row, UAV_allometry_join,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
#this may be helpful. The small block is represented as branchier
UAV_allometry_join$branchiness <- branchiness_lm_2022$residuals

write.csv(UAV_allometry_join, file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_allometry_processed.csv")

library(ggmap)
library(terra)
library(sf)
nc <- read_sf(system.file("gpkg/nc.gpkg", package = "sf"))
plot(nc["BIR74"], reset = FALSE, key.pos = 4)
plot(st_buffer(nc[1,1], units::set_units(10, km)), col = 'NA', 
     border = 'red', lwd = 2, add = TRUE)

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

#spatial with GIS dat
UAV_allometry_sf <- st_as_sf(UAV_allometry_join, wkt = "geometry")

str(UAV_join)
ggplot() + 
  geom_sf(data = UAV_allometry_sf, aes(fill = branchiness))+
  scale_y_continuous(breaks = -1:1)+
  scale_fill_gradient()
  

