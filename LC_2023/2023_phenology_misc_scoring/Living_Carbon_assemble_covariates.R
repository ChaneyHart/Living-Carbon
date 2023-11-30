#2023_field_season_covariates

#read in growth dat

growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
summary(growth$block)
growth %>% group_by(block) %>% tally()

SPAD_dat <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/SPAD_full.csv")

UAV_dat_2022 <- read.csv(file = "LC_2022/2022_UAV_tables/LC_full_complete.csv")

Leaf_morphology_dat_2022 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/2022_leaf_morphology/SLA_2022_raw.csv")

covariates_list <- list(SPAD_dat,UAV_dat_2022,Leaf_morphology_dat_2022)

covariates_inner_join <- covariates_list %>% reduce(inner_join, by = "ID")

covariates_full_join <- covariates_list %>% reduce(full_join, by = "ID")

write.csv(covariates_inner_join, file = "LC_2023/2023_growth_inventory_analysis/covariates_inner.csv")
write.csv(covariates_full_join, file = "LC_2023/2023_growth_inventory_analysis/covariates_full.csv")

