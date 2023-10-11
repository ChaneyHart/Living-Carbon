library(dplyr)

#read in most up to date inventory
#read in most up to date inventory

LC_inventory <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
aug_7_inv <- read.csv(file = "LC_2023/8_7_inventory.csv")
aug_7_inv <- subset(aug_7_inv, select = c(ID,Drought.Score))
LC_inventory <- inner_join(LC_inventory, aug_7_inv, by = "ID")


LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))

#Refine population from which to sample from
#Want to select from events

LC_diurnal_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")
#confine to large block where trees have more lateral branches

#Do not consider stunted trees - trees whose height is 1 stdev below mean
mean(LC_diurnal_pop$H497)
sd(LC_diurnal_pop$H497)
hist(LC_diurnal_pop$H497)
ht_threshold <- mean(LC_diurnal_pop$H497) - (1*(sd(LC_diurnal_pop$H497)))
LC_diurnal_pop <- subset(LC_diurnal_pop, H497 >= ht_threshold)

#remove trees in border
LC_diurnal_pop <-subset(LC_diurnal_pop, row != 1)
LC_diurnal_pop <-subset(LC_diurnal_pop, column != 30)

#Remove known damaged trees
LC_diurnal_pop <-subset(LC_diurnal_pop, ID != "LCOR-197")

##
#exclude irrigation set
Irrigation_control_set <- read.csv("LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
LC_diurnal_pop <- LC_diurnal_pop[!(LC_diurnal_pop$ID %in% Irrigation_set$ID),]


#sample droughted trees
LC_diurnal_sample_total_extras <- LC_diurnal_pop %>% group_by(event_short) %>% sample_n(12)
LC_diurnal_sample_total <- LC_diurnal_pop %>% group_by(event_short) %>% sample_n(10)
LC_diurnal_sample_total_extras$method <- "Li600/MiniPAM"
LC_diurnal_sample_total_extras$treatment <- "drought"

LC_diurnal_sample_6800 <- LC_diurnal_sample_total %>% group_by(event_short) %>% sample_n(3)
LC_diurnal_sample_6800$method <- "Li6800"
LC_diurnal_sample_6800$treatment <- "drought"

#sample control trees
LC_diurnal_control <- subset(Irrigation_control_set, select = c(row, column,ID,event,event_short,block,H497))
LC_diurnal_control$method <- "Li-600/MiniPAM"
LC_diurnal_control$treatment <- "control"
LC_diurnal_control_6800 <- LC_diurnal_control %>% group_by(event_short) %>% sample_n(2)
LC_diurnal_control_6800$method <- "Li6800"
LC_diurnal_control_6800$treatment <- "control"

LC_diurnal_drought_sample <- rbind(LC_diurnal_sample_total_extras,LC_diurnal_sample_6800)
LC_diurnal_control_sample <- rbind(LC_diurnal_control, LC_diurnal_control_6800)

LC_diurnal_sample <- rbind(LC_diurnal_drought_sample, LC_diurnal_control_sample)
#write.csv(LC_diurnal_sample_6800, file = "LC_2023/LC_diurnal_Li6800_6_14.csv")
#write.csv(LC_diurnal_sample_total, file = "LC_2023/LC_diurnal_total_6_14.csv")

#write.csv(LC_diurnal_sample_6800, file = "LC_2023/LC_diurnal_Li6800_6_26.csv")
#write.csv(LC_diurnal_sample_total, file = "LC_2023/LC_diurnal_total_6_26.csv")

write.csv(LC_diurnal_sample_6800, file = "LC_2023/LC_diurnal_Li6800_7_20.csv")
write.csv(LC_diurnal_sample_total, file = "LC_2023/LC_diurnal_total_7_20.csv")
write.csv(LC_diurnal_sample_total_extras, file = "LC_2023/LC_diurnal_total_extras_7_20.csv")

write.csv(LC_diurnal_sample, file = "LC_diurnal_sample_8_29.csv")



#Sample for diurnal analysis w/ non-photorespiratory conditions

LC_lowO2_drought <- LC_diurnal_pop %>% group_by(event_short) %>% sample_n(1)
LC_lowO2_control <- LC_diurnal_control %>% group_by(event_short) %>% sample_n(1)

LC_lowO2 <- rbind(LC_lowO2_drought,LC_lowO2_control)

write.csv(LC_lowO2, file = "LC_lowO2_sample.csv")

#Sampling for drone flight

LC_drone_sample_total <- LC_diurnal_pop %>% group_by(event_short) %>% sample_n(8)
LC_drone_sample <- LC_drone_sample_total %>% group_by(event_short) %>% sample_n(5)

write.csv(LC_drone_sample, file = "LC_2023/LC_drone_sample.csv")
write.csv(LC_drone_sample_total, file = "LC_2023/LC_drone_sample_total.csv")
