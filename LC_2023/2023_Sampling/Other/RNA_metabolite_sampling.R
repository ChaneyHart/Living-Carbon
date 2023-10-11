# selecting random trees for RNA sampling
library(purrr)
library(tidyr)
library(dplyr)

#read in most up to date inventory

LC_inventory <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
aug_7_inv <- read.csv(file = "LC_2023/8_7_inventory.csv")
aug_7_inv <- subset(aug_7_inv, select = c(ID,Drought.Score))
LC_inventory <- inner_join(LC_inventory, aug_7_inv, by = "ID")

LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))
#add info about section of field


#Refine population from which to sample from
#Want to select from events
LC_RNA_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_RNA_pop$H497)
sd(LC_RNA_pop$H497)
ht_threshold <- mean(LC_RNA_pop$H497) - (1*(sd(LC_RNA_pop$H497)))
LC_RNA_pop <- subset(LC_RNA_pop, H497 >= ht_threshold)

LC_RNA_pop <-subset(LC_RNA_pop, row != 1)
LC_RNA_pop <-subset(LC_RNA_pop, column != 30)


#how many of each event are there?
LC_RNA_test <- LC_RNA_pop %>% group_by(event_short)
group_data(LC_RNA_test)

#Two independent samples for early and late season sampling
#If necessary, in the end - can pool
RNA_pop_1 <- LC_RNA_pop %>% group_by(event_short) %>% sample_n(8)

RNA_pop_11 <- LC_RNA_pop %>% group_by(event_short) %>% sample_n(1)

RNA_pop_2 <- LC_RNA_pop %>% group_by(event_short) %>% sample_n(8)

RNA_sample_1 <- RNA_pop_1 %>% group_by(event_short) %>% sample_n(6)

RNA_sample_2 <- RNA_pop_2 %>% group_by(event_short) %>% sample_n(6)

write.csv(RNA_sample_1, file = "RNA_sample_early_2023.csv")

#write.csv(RNA_sample_2, file = "RNA_sample_late_2023.csv")

## Sampling metabolites from the population of trees being sampled for RNA expression


metabolite_sample_1 <- subset(RNA_pop_1, event_short == "13-15E" | event_short == "2H" | event_short == "5A" |  event_short == "16-20" | event_short == "8-9D" )
# select other two samples for total of 10 per event

LC_metabolite_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")
#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_metabolite_pop$H497)
sd(LC_metabolite_pop$H497)
ht_threshold <- mean(LC_metabolite_pop$H497) - (1*(sd(LC_metabolite_pop$H497)))
LC_metabolite_pop <- subset(LC_metabolite_pop, H497 >= ht_threshold)

LC_metabolite_pop <-subset(LC_metabolite_pop, row != 1)
LC_metabolite_pop <-subset(LC_metabolite_pop, column != 30)
LC_metabolite_pop <- subset(LC_metabolite_pop, event_short == "13-15E" | event_short == "2H" | event_short == "5A" |  event_short == "16-20" | event_short == "8-9D")

metabolite_sample_1.1 <- LC_metabolite_pop %>% group_by(event_short) %>% sample_n(3)
metabolite_sample_1_final <- unique(rbind(metabolite_sample_1, metabolite_sample_1.1))

#metabolite_sample_2 <- subset(RNA_pop_2, event_short == "13-15E" | event_short == "2H" | event_short == "5A" |  event_short == "16-20" | event_short == "8-9D" )



write.csv(metabolite_sample_1_final, file = "Metabolite_sample_early_2023.csv")

#write.csv(metabolite_sample_2, file = "Metabolite_sample_late_2023.csv")



###Metabolite sampling
LC_metabolite_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_metabolite_pop$H497)
sd(LC_metabolite_pop$H497)
ht_threshold <- mean(LC_metabolite_pop$H497) - (1*(sd(LC_metabolite_pop$H497)))
LC_metabolite_pop <- subset(LC_metabolite_pop, H497 >= ht_threshold)

LC_metabolite_pop <-subset(LC_metabolite_pop, row != 1)
LC_metabolite_pop <-subset(LC_metabolite_pop, column != 30)
LC_metabolite_pop <- subset(LC_metabolite_pop, event_short == "13-15E" | event_short == "2H" | event_short == "5A" |  event_short == "16-20" | event_short == "8-9D")

#want metabolite samples to overlap with those taken for RNA + 4 more trees



metabolite_sample_1 <- LC_metabolite_pop %>% group_by(event_short) %>% sample_n(10)

metabolite_sample_2 <- LC_metabolite_pop %>% group_by(event_short) %>% sample_n(10)

write.csv(metabolite_sample_1, file = "Metabolite_sample_early_2023.csv")

write.csv(metabolite_sample_2, file = "Metabolite_sample_late_2023.csv")


metabolite_sample_corrected <- read.csv("Metabolite_sample_early_2023.csv")

metabolite_early_1 <- subset(metabolite_sample_corrected, backup != 'yes') %>% group_by(event_short) %>% sample_n(5)

write.csv(metabolite_early_1, file = "metabolite_sample_7_12_23.csv")



## re-picking sample for late season sampling

LC_inventory <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
aug_7_inv <- read.csv(file = "LC_2023/8_7_inventory.csv")
aug_7_inv <- subset(aug_7_inv, select = c(ID,Drought.Score))
LC_inventory <- inner_join(LC_inventory, aug_7_inv, by = "ID")

LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497","Drought.Score"))
#add info about section of field


#Refine population from which to sample from
#Want to select from events
LC_RNA_pop_late <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_RNA_pop_late$H497)
sd(LC_RNA_pop_late$H497)
ht_threshold <- mean(LC_RNA_pop_late$H497) - (1*(sd(LC_RNA_pop_late$H497)))
LC_RNA_pop <- subset(LC_RNA_pop_late, H497 >= ht_threshold)

LC_RNA_pop_late <-subset(LC_RNA_pop_late, row != 1)
LC_RNA_pop_late <-subset(LC_RNA_pop_late, column != 30)



#exclude samples in early sampling

early_RNA <- read.csv(file = "LC_2023/2023_Sampling/Other/RNA_sample_early_2023.csv")
LC_RNA_pop_late <- LC_RNA_pop_late[!(LC_RNA_pop_late$ID %in% early_RNA$ID),]

early_metabolite<- read.csv(file = "LC_2023/2023_Sampling/Other/Metabolite_sample_early_2023.csv")
LC_RNA_pop_late <- LC_RNA_pop_late[!(LC_RNA_pop_late$ID %in% early_metabolite$ID),]

#exclude irrigation set
Irrigation_set <- read.csv("LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
LC_RNA_pop_late <- LC_RNA_pop_late[!(LC_RNA_pop_late$ID %in% Irrigation_set$ID),]

#Do not consider trees with drought score 3 or higher
str(LC_RNA_pop_late)
LC_RNA_pop_late$Drought.Score <- as.numeric(LC_RNA_pop_late$Drought.Score)
LC_RNA_pop_late <- subset(LC_RNA_pop_late, Drought.Score < 3)


#how many of each event are there?
LC_RNA_test <- LC_RNA_pop_late %>% group_by(event_short)
group_data(LC_RNA_test)

#narrowing down the ones that can be picked from the early sample
early_RNA2 <- early_RNA[!(early_RNA$ID %in% Irrigation_set$ID),]
Lc_drought_score <- subset(LC_sample_list, select = c(ID,Drought.Score))
early_RNA2 <- inner_join(early_RNA2, Lc_drought_score, by = "ID")
early_RNA2 <- subset(early_RNA2, Drought.Score < 3)

#Final late sample will be three from previous sample and 3 from new sample

Late_RNA_sample_a <- LC_RNA_pop_late %>% group_by(event_short) %>% sample_n(3)
Late_RNA_sample_a <- subset(Late_RNA_sample_a, select = c(row,column,ID,event_short,H497,Drought.Score))
Late_RNA_sample_b <- early_RNA2 %>% group_by(event_short) %>% sample_n(3)
Late_RNA_sample_b <- subset(Late_RNA_sample_b, select = c(row,column,ID,event_short,H497,Drought.Score))
Late_RNA_sample_b$Drought.Score <- as.double(Late_RNA_sample_b$Drought.Score)

Late_RNA_sample <- rbind(Late_RNA_sample_a,Late_RNA_sample_b)

write.csv(Late_RNA_sample, file = "LC_2023/2023_Sampling/Other/RNA_sample_final_late_8_14_23.csv")

##Metabolite sampling

#the list of ones from the early sample that are shared between RNA and metabolites and expected to still be viable
early_metabolite2 <- subset(early_metabolite, select = ID)

early_RNA3 <- early_RNA[!(early_RNA$ID %in% Irrigation_set$ID),]

early_RNA3 <- inner_join(early_RNA3, Lc_drought_score, by = "ID")
early_RNA3 <- subset(early_RNA3, Drought.Score < 4)

early_sampling_combined <- inner_join(early_RNA3, early_metabolite2)


#narrowing down the ones that can be picked from the early sample
early_metabolite3 <- early_metabolite[!(early_metabolite$ID %in% early_sampling_combined$ID),]
early_metabolite3 <- early_metabolite3[!(early_metabolite3$ID %in% Irrigation_set$ID),]
early_metabolite3 <- inner_join(early_metabolite3, Lc_drought_score, by = "ID")
early_metabolite3 <- subset(early_metabolite3, Drought.Score < 4)

LC_H497 <- subset(LC_sample_list, select = c(ID,H497))
early_metabolite3 <- inner_join(early_metabolite3, LC_H497)

#the list of trees not sampled in early sampling
late_metabolite_pop <-  subset(LC_RNA_pop_late, event_short == "13-15E" | event_short == "2H" | event_short == "5A" |  event_short == "16-20" | event_short == "8-9D" )

#sampling will be from three populations:
#trees that were sampled for RNA and metabolites in July (as many as possible)
#Other trees sampled for metabolites but not RNA in July (as needed)
#Filling out to sample size of 10 with trees not sampled in July
Late_metabolite_sample_a <- early_sampling_combined %>% group_by(event_short) %>% sample_n(4)
Late_metabolite_sample_a <- subset(Late_metabolite_sample_a, select = c(row,column,ID,event_short,H497,Drought.Score))
Late_metabolite_sample_a$Drought.Score <- as.double(Late_metabolite_sample_a$Drought.Score)
Late_metabolite_sample_b <- early_metabolite3 %>% group_by(event_short) %>% sample_n(4)
Late_metabolite_sample_b <- subset(Late_metabolite_sample_b, select = c(row,column,ID,event_short,H497,Drought.Score))
Late_metabolite_sample_b$Drought.Score <- as.double(Late_metabolite_sample_b$Drought.Score)
Late_metabolite_sample_c <- late_metabolite_pop %>% group_by(event_short) %>% sample_n(2)
Late_metabolite_sample_c <- subset(Late_metabolite_sample_c, select = c(row,column,ID,event_short,H497,Drought.Score))                                                             
Late_metabolite_sample_c$Drought.Score <- as.double(Late_metabolite_sample_c$Drought.Score)

Late_metabolite_sample <- rbind(Late_metabolite_sample_a,Late_metabolite_sample_b,Late_metabolite_sample_c)

write.csv(Late_metabolite_sample, file = "LC_2023/2023_Sampling/Metabolite_sample_final_late_8_16_23.csv")
    
Late_metabolite_sample

Late_metabolite_WP_sample <- Late_metabolite_sample %>% group_by(event_short) %>% sample_n(5)                                                         
