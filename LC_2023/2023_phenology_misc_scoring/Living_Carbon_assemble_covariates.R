#2023_field_season_covariates
library(dplyr)
library(ggplot2)
library(sf)
library(tidyverse)

#read in growth dat

growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")

# read in covariates

SPAD_dat <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/SPAD_summary.csv")

UAV_indices <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_indices.csv")
UAV_allometry <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_allometry_processed.csv")

#Leaf_morphology_dat_2022 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/2022_leaf_morphology/SLA_2022_raw.csv")

covariates_list <- list(SPAD_dat,UAV_indices,UAV_allometry)

covariates_inner_join <- covariates_list %>% reduce(inner_join, by = "ID")

covariates_full_join <- covariates_list %>% reduce(full_join, by = "ID")

write.csv(covariates_inner_join, file = "LC_2023/2023_growth_inventory_analysis/covariates_inner.csv")
write.csv(covariates_full_join, file = "LC_2023/2023_growth_inventory_analysis/covariates_full.csv")



growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")

#only intersted in trees matching growth dat
covariates_set <- covariates_inner_join[(covariates_inner_join$ID %in% growth$ID),]
# remove covariates with incomplete coverage

#remove extraneous variables
covariates_set_lite <- subset(covariates_set, select = -c(X.x,row.x,column.x,X.y,row.y,column.y,Construct.x,Event.x,area.x,geometry.x,X.1,X,Construct.y,Event.y,area.y,geometry.y,block.x,block.y,row,column,shadow_pix,blue_mean,blue_std,green_mean,green_std,red_mean,red_std,nir_mean,nir_std,re_mean,re_std,TGI_std,GRVI_std,MGRVI_std,EXG_std,EXGR_std,NDVI_std,NDRE_std,GNDVI_std,SAVI_std,OSAVI_std,MSAVI_std,GCI_std,RECI_std,event_short))

#look at correlations

write.csv(covariates_set_lite,file = "LC_2023/2023_growth_inventory_analysis/processed/covariates_set.csv")            

covariates_cor_set<- na.omit(covariates_set_lite)
covariates_cor <- cor((covariates_cor_set[,-1]), use = "pairwise.complete.obs", method = "pearson")

png(filename = "LC_2023/2023_growth_inventory_analysis/figures/full_covariate_cor.png",height = 8, width =8, units = "in",res = 300)
corrplot(as.matrix(covariates_cor[]))
dev.off()

#SPAD

library(latticeExtra)

png(filename = "LC_2023/2023_growth_inventory_analysis/figures/SPAD_june_2022.png",width = 8,height=4, units = "in",res = 300)
levelplot(covariates_set$SPAD_June_2022 ~-column*row, covariates_set,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=1000))
dev.off()

png(filename = "LC_2023/2023_growth_inventory_analysis/figures/vol_spatial.png",width = 8,height=4, units = "in",res = 300)
levelplot(covariates_set$V801 ~-column*row, covariates_set,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=1000))
dev.off()

png(filename = "LC_2023/2023_growth_inventory_analysis/figures/Aug_spad_2022.png",width = 8,height=4, units = "in",res = 300)
levelplot(covariates_set$SPAD_aug_2022 ~-column*row, covariates_set,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=1000))
dev.off()

SPAD_vol_scatter <- ggplot(covariates_set,aes(x=SPAD_June_2022,y=V801))+
  geom_point()+
  ylab("Volume index (cm3) in Sep 2023")+
  xlab("Chl density (SPAD) June 2022")+
  theme_bw()
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/figures/SPAD_vol_scatter.png",SPAD_vol_scatter, dpi=300)




