#Sampling for REsponse curves


library(dplyr)
library(readxl)

#read in most up to date inventory

LC_inventory <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
#read in Drought_scores
aug_7_inv <- read.csv(file = "LC_2023/8_7_inventory.csv")
aug_7_inv <- subset(aug_7_inv, select = c(ID,Drought.Score))
LC_inventory <- inner_join(LC_inventory, aug_7_inv, by = "ID")
LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497","Drought.Score"))

#Refine population from which to sample from
#Want to select from events

LC_phys_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_phys_pop$H497)
sd(LC_phys_pop$H497)
ht_threshold <- mean(LC_phys_pop$H497) - (1.5*(sd(LC_phys_pop$H497)))
LC_phys_pop_II <- subset(LC_phys_pop, H497 >= ht_threshold)

#Do not consider trees already sampled

Sampled_trees <- read_excel(path = "LC_2023/Phys_running_list.xlsx")

LC_phys_pop_III <- LC_phys_pop_II[!(LC_phys_pop_II$ID %in% Sampled_trees$ID),]

#do not consider water bag trees

Irrigation_set <- read.csv("LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
LC_phys_pop_III <- LC_phys_pop_III[!(LC_phys_pop_III$ID %in% Irrigation_set$ID),]

#Do not consider trees with drought score 3 or higher
str(LC_phys_pop_III)
LC_phys_pop_III$Drought.Score <- as.numeric(LC_phys_pop_III$Drought.Score)
LC_phys_pop_VI <- subset(LC_phys_pop_III, Drought.Score < 3)

LC_phys_pop_V <- LC_phys_pop_VI %>% group_by(event_short)
group_data(LC_phys_pop_V)

LC_phys_sample <- LC_phys_pop_VI%>% group_by(event_short) %>% sample_n(3)


# randomly sample these so that the order of trees is not alphanumeric
#LC_phys_sample_1 <- LC_phys_sample %>% group_by(event_short) %>% sample_n(1, replace = TRUE)

#LC_phys_sample_2 <- LC_phys_sample %>% group_by(event_short) %>% sample_n(1, replace = TRUE)

#LC_phys_sample_3 <- LC_phys_sample %>% group_by(event_short) %>% sample_n(1, replace = TRUE)

#write.csv(LC_phys_sample, file = "LC_2023/phys_sample_6_14.csv")
#write.csv(LC_phys_sample, file = "LC_2023/phys_sample_6_22.csv")
#write.csv(LC_phys_sample, file = "LC_2023/phys_sample_7_20.csv")
#write.csv(LC_phys_sample, file = "LC_2023/phys_sample_7_25.csv")
#write.csv(LC_phys_sample, file = "LC_2023/phys_sample_8_9.csv")
write.csv(LC_phys_sample, file = "LC_2023/phys_sample_8_25.csv")
