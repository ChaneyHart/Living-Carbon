# selecting random trees to take cuttings from for living carbon greenhouse planting.


library(dplyr)

#read in most up to date inventory

LC_inventory <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))

#Refine population from which to sample from
#Want to select from events

LC_cutting_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")
#confine to large block where trees have more lateral branches
LC_cutting_pop <- subset(LC_cutting_pop, block == "large")
#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_cutting_pop$H497)
sd(LC_cutting_pop$H497)
ht_threshold <- mean(LC_cutting_pop$H497) - (2*(sd(LC_cutting_pop$H497)))
LC_cutting_pop_II <- subset(LC_cutting_pop, H497 >= ht_threshold)

set.seed(68)
LC_cutting_sample <- LC_cutting_pop_II%>% group_by(event_short) %>% sample_n(8)



write.csv(LC_cutting_sample, file = "LC_2023/cutting_sample.csv")



