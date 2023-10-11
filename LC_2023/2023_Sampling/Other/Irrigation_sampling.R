# selecting random trees to irrigate from for living carbon planting.

library(purrr)
library(tidyr)
library(dplyr)

#read in most up to date inventory

LC_inventory <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))

#add info about section of field

LC_sample_list <- LC_sample_list %>% mutate(section = case_when(
  column <= 8 & row <= 8 ~ 1,
  column <= 8 & row > 8 ~ 2,
  column > 8 & column <= 16 & row <= 8 ~ 3,
  column > 8 & column <= 16 & row > 8 ~ 4,
  column > 16 & column <= 24 & row <= 8 ~5,
  column > 16 & column <= 24 & row > 8 ~ 6,
  column > 24 & column <= 28 & row <= 8 ~7,
  column > 24 & column <= 28 & row > 8 ~8,
  column > 44 & column <= 50 & row <=8 ~9,
  column > 44 & column <= 50 & row > 8 ~10,
  column > 50 & column < 56 & row <= 8 ~ 11,
  column > 50 & column < 56 & row > 8 ~ 12
))

#Refine population from which to sample from
#Want to select from events
LC_irrigation_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#Do not consider stunted trees - trees whose height is 2 stdev below mean
mean(LC_irrigation_pop$H497)
sd(LC_irrigation_pop$H497)
ht_threshold <- mean(LC_irrigation_pop$H497) - (1*(sd(LC_irrigation_pop$H497)))
LC_irrigation_pop <- subset(LC_irrigation_pop, H497 >= ht_threshold)



#confine to accessible part of field


LC_irrigation_pop <- subset(LC_irrigation_pop, section == "2" | section == "4")



LC_irrigation_sample <- LC_irrigation_pop %>% group_by(event_short) %>% sample_n(3)


write.csv(LC_irrigation_sample, file = "LC_irr_sample_large.csv")
## a solution with unequal sample sizes

iris %>%
  group_by(Species) %>% 
  nest() %>%            
  ungroup() %>% 
  mutate(n = c(2, 5, 3)) %>% 
  mutate(samp = map2(data, n, sample_n)) %>% 
  select(-data) %>%
  unnest(samp)


LC_irrigation_pop_test <- LC_irrigation_pop %>% group_by(event_short)
group_data(LC_irrigation_pop_test)


LC_irrigation_popII <- LC_irrigation_pop %>%
  group_by(event_short) %>% 
  nest() %>%            
  ungroup() %>% 
  mutate(n = c(2, 3, 2, 2, 3, 3, 2, 2)) %>% 
  mutate(samp = map2(data, n, sample_n)) %>% 
  select(-data) %>%
  unnest(samp)

write.csv(LC_irrigation_popII, file="Lc_irrigation_sampleII.csv")


LC_inventory_test <- LC_inventory %>% group_by(event_short)
group_data(LC_inventory_test)
