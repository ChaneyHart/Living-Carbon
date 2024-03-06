#data wrangling

library(dplyr)
library(ggplot2)


#read in growth dat

growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")
growth <- subset(growth, select = c(X,ID,row,column,event,event_short,construct,construct2,block,H49,H801,D49,D801,V49,V801))

# read in covariates if they rep trees in growth set

SPAD_dat <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/SPAD_summary.csv")
SPAD_dat <- SPAD_dat %>% filter(ID %in% growth$ID)
UAV_indices <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_indices.csv")
UAV_indices <- subset(UAV_indices, select = -c(Construct,Event,area,blue_mean,blue_std,red_mean,red_std,green_mean,green_std,nir_std,re_std,TGI_std,GRVI_std,MGRVI_std,EXG_std,EXGR_std,NDVI_std,NDRE_std,GNDVI_std,SAVI_std,OSAVI_std,MSAVI_std,GCI_std,RECI_std))
UAV_indices <- UAV_indices %>% filter(ID %in% growth$ID)
UAV_allometry <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/UAV_allometry_processed.csv")
UAV_allometry <- subset(UAV_allometry, select = c(ID,shadow_area,branchiness))

UAV_allometry <- UAV_allometry %>% filter(ID %in% growth$ID)

Leaf_area <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/2024_leaf_mass_area.csv")
Leaf_area <- subset(Leaf_area, select = c(ID,mass_avg,area_avg,SLA))
Leaf_area <- Leaf_area %>% filter(ID %in% growth$ID)

#these datasets do not all match in size. For now, the sample of the trees that were sampled for leaf area will be analyzed

growth_reduced <- growth %>% filter(ID %in% Leaf_area$ID)
SPAD_reduced <- SPAD_dat %>% filter(ID %in% Leaf_area$ID)
UAV_indices_reduced <- UAV_indices %>% filter(ID %in% Leaf_area$ID)
UAV_allometry_reduced <- UAV_allometry %>% filter(ID %in% Leaf_area$ID)

#merge datasets
str(growth_reduced)
str(SPAD_reduced)
str(Leaf_area)
str(UAV_indices_reduced)

full_dat <- inner_join(growth_reduced,SPAD_reduced, by = c("ID","row","column"))
full_dat <- inner_join(full_dat,UAV_indices_reduced, by = c("ID","row","column"))
full_dat <- inner_join(full_dat,UAV_allometry_reduced, by = c("ID"))
full_dat <- inner_join(full_dat,Leaf_area, by = c("ID"))

full_dat <- full_dat %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

processed_growth_dat <- subset(full_dat, select = -c(X.x,X.y,X,event,construct,construct2))


write.csv(processed_growth_dat,file = "LC_2023/2023_growth_inventory_analysis/processed/processed_growth_dat.csv")

