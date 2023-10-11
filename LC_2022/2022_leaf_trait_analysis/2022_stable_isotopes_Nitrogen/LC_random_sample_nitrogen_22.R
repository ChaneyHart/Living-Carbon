# selecting random trees for living carbon analyses


library(dplyr)

#read in most up to date inventory

LC_inventory <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block"))

#Refine population from which to sample from
LC_Nitrogen_pop <- subset(LC_sample_list, block == "large")

LC_Nitrogen_pop <- subset(LC_Nitrogen_pop, event_short != "4B")
LC_Nitrogen_pop <- subset(LC_Nitrogen_pop, event_short != "1")
#check correct events are present
count(LC_Nitrogen_pop$event)

set.seed(55)
LC_Nitrogen_sample <- LC_Nitrogen_pop%>% group_by(event_short) %>% sample_n(11)



write.csv(LC_Nitrogen_sample, file = "2022_leaf_trait_analysis/2022_stable_isotopes_Nitrogen/LC_Nitrogen_isotope_datasheet_2022.csv")
