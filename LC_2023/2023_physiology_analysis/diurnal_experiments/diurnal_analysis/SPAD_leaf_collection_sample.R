library(dplyr)

#read in most up to date inventory
#read in most up to date inventory

LC_inventory <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
aug_7_inv <- read.csv(file = "LC_2023/8_7_inventory.csv")
aug_7_inv <- subset(aug_7_inv, select = c(ID,Drought.Score))
LC_inventory <- inner_join(LC_inventory, aug_7_inv, by = "ID")


Leaf_collection_datasheet <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))

Leaf_collection_datasheet$Leaf_collection <- ifelse(LC_sample_list$event_short == "13-15E"|LC_sample_list$event_short == "5C"|LC_sample_list$event_short == "5A"|LC_sample_list$event_short == "4A"|LC_sample_list$event_short == "2H"|LC_sample_list$event_short == "16-20"|LC_sample_list$event_short == "8-9D"|LC_sample_list$event_short == "CT3", "yes","no")

write.csv(Leaf_collection_datasheet, file = "LC_2023/2023_growth_inventory_analysis/leaf_collection_scoring.csv")
















