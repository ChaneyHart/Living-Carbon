# selecting trees for pre-dawn water potential sampling from for living carbon greenhouse planting.


library(dplyr)

#read in most up to date inventory

LC_inventory <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_sample_list <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","H497"))

#define blocks 



LC_WP_pop <- LC_sample_list %>% mutate(section = case_when(
  column <= 8 & row <= 8 ~ 1,
  column <= 8 & row > 8 ~ 2,
  column > 8 & column <= 16 & row <= 8 ~ 3,
  column > 8 & column <= 16 & row > 8 ~ 4,
  column > 16 & column <= 24 & row <= 8 ~5,
  column > 16 & column <= 24 & row > 8 ~ 6,
  column > 24 & column <= 30 & row <= 8 ~7,
  column > 24 & column <= 30 & row > 8 ~8,
  column > 44 & column <= 50 & row <=8 ~9,
  column > 44 & column <= 50 & row > 8 ~10,
  column > 50 & column < 56 & row <= 8 ~ 11,
  column > 50 & column < 56 & row > 8 ~ 12
))


mean(LC_WP_pop$H497)
sd(LC_WP_pop$H497)
ht_threshold <- mean(LC_WP_pop$H497) - (2*(sd(LC_WP_pop$H497)))
LC_WP_popII <- subset(LC_WP_pop, H497 >= ht_threshold)

LC_WP_sample <- LC_WP_popII %>% group_by(section) %>% sample_n(3)
write.csv(LC_WP_sample, file = "8_2_23_predawnWP_datasheet.csv")


LC_fv_fm_pop <- subset(LC_sample_list, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

mean(LC_fv_fm_pop$H497)
sd(LC_fv_fm_pop$H497)
ht_threshold <- mean(LC_fv_fm_pop$H497) - (1*(sd(LC_fv_fm_pop$H497)))
LC_fv_fm_popII <- subset(LC_fv_fm_pop, H497 >= ht_threshold)

LC_fv_fm_popII <- LC_fv_fm_popII %>% group_by(event_short)
group_data(LC_fv_fm_popII)

LC_Fv_fm_sample <- LC_fv_fm_pop %>% group_by(event_short) %>% sample_n(19)

write.csv(LC_Fv_fm_sample, file = "8_1_23_Fv_Fm_datasheet.csv")
