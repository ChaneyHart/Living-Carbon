library(ggplot2)
library(dplyr)
library(readxl)

#read in data

growth_dat_8_7 <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/8_7_inventory.csv")

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("ID","event","event_short","block","construct","construct2","H497"))

LC_meta <- LC_meta %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

growth_dat_8_7 <- inner_join(LC_meta, growth_dat_8_7, by = "ID")


growth_dat_8_7 <- growth_dat_8_7 %>% mutate(section = case_when(
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

growth_dat_8_7$section <- as.factor(growth_dat_8_7$section)

growth_dat_8_7$d0 <- ifelse(growth_dat_8_7$Drought.Score == "0", 1, 0)
growth_dat_8_7$d1 <- ifelse(growth_dat_8_7$Drought.Score == "1", 1, 0)
growth_dat_8_7$d2 <- ifelse(growth_dat_8_7$Drought.Score == "2", 1, 0)
growth_dat_8_7$d3 <- ifelse(growth_dat_8_7$Drought.Score == "3", 1, 0)
growth_dat_8_7$d4 <- ifelse(growth_dat_8_7$Drought.Score == "4", 1, 0)
growth_dat_8_7$d5 <- ifelse(growth_dat_8_7$Drought.Score == "5", 1, 0)

str(growth_dat)
drought_event_8_7 <- growth_dat %>% group_by(event_short) %>% summarise(
  d0 = sum(d0),
  d1 = sum(d1),
  d2 = sum(d2),
  d3 = sum(d3),
  d4 = sum(d4),
  d5 = sum(d4)
)

drought_event_8_7$date <- "Aug_7"

drought_section_8_7 <- growth_dat %>% group_by(section) %>% summarise(
  d0 = sum(d0),
  d1 = sum(d1),
  d2 = sum(d2),
  d3 = sum(d3),
  d4 = sum(d4),
  d5 = sum(d4)
)

drought_section_8_7$date <- "Aug_7"


#Drought scoring end of Aug

drought_8_21 <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/8_21_drought_score_SPAD.csv")


LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("ID","event","event_short","block","construct","construct2","H497"))

LC_meta <- LC_meta %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

growth_dat_8_21 <- inner_join(LC_meta, drought_8_21, by = "ID")


growth_dat_8_21 <- growth_dat_8_21 %>% mutate(section = case_when(
  column <= 8 & row <= 8 ~ 1,
  column <= 8 & row > 8 ~ 2,
  column > 8 & column <= 16 & row <= 8 ~ 3,
  column > 8 & column <= 16 & row > 8 ~ 4,
  column > 16 & column <= 24 & row <= 8 ~5,
  column > 16 & column <= 24 & row > 8 ~ 6,
  column > 24 & column <= 28 & row <= 8 ~ 7,
  column > 24 & column <= 28 & row > 8 ~8,
  column > 44 & column <= 50 & row <=8 ~9,
  column > 44 & column <= 50 & row > 8 ~10,
  column > 50 & column < 56 & row <= 8 ~ 11,
  column > 50 & column < 56 & row > 8 ~ 12
))

growth_dat_8_21$section <- as.factor(growth_dat_8_21$section)

growth_dat_8_21$d0 <- ifelse(growth_dat_8_21$Drought.Score == "0", 1, 0)
growth_dat_8_21$d1 <- ifelse(growth_dat_8_21$Drought.Score == "1", 1, 0)
growth_dat_8_21$d2 <- ifelse(growth_dat_8_21$Drought.Score == "2", 1, 0)
growth_dat_8_21$d3 <- ifelse(growth_dat_8_21$Drought.Score == "3", 1, 0)
growth_dat_8_21$d4 <- ifelse(growth_dat_8_21$Drought.Score == "4", 1, 0)
growth_dat_8_21$d5 <- ifelse(growth_dat_8_21$Drought.Score == "5", 1, 0)

str(growth_dat)
drought_event_8_21 <- growth_dat_8_21 %>% group_by(event_short) %>% summarise(
  d0 = sum(d0),
  d1 = sum(d1),
  d2 = sum(d2),
  d3 = sum(d3),
  d4 = sum(d4),
  d5 = sum(d4)
)

drought_section_8_21 <- growth_dat_8_21 %>% group_by(section) %>% summarise(
  d0 = sum(d0),
  d1 = sum(d1),
  d2 = sum(d2),
  d3 = sum(d3),
  d4 = sum(d4),
  d5 = sum(d4)
)

drought_event_8_21$date <- "Aug_21"
drought_section_8_21$date <- "Aug_21"

drought_combined <- rbind(drought_event_8_7, drought_event_8_21)

str(drought_combined)
ggplot(drought_combined, aes(x=date, color = event_short))+
  geom_point(aes(y = d0))+
  geom_line(aes(y = d0))+
  geom_point(aes(y = d1))+
  geom_line(aes(y = d1))+
  geom_point(aes(y = d2))+
  geom_point(aes(y = d3))+
  geom_point(aes(y = d4))+
  geom_point(aes(y = d5))

