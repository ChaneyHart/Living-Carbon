library(ggplot2)
library(dplyr)
library(readxl)

#read in data

growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/05_06_DBHH.csv", skip = 5)
growth_dat$Vol_index <- 0.333*((3.14*((growth_dat$DBH/2)^2))*growth_dat$Height)

excluded_trees_raw <- read_xlsx("LC_2022/2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/LC_all trees to exclude.xlsx",1)
excluded_trees <- as.vector(excluded_trees_raw$`to exclude - dead, replaced, small, diseased`)

biomass_set <- filter(growth_dat, !ID %in% excluded_trees)
biomass_set <- filter(biomass_set, DBH > 0)


LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2"))

LC_meta <- LC_meta %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

LC_meta <- LC_meta %>% mutate(section = case_when(
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


LC_meta <- LC_meta %>% mutate(section = case_when(
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


biomass_set <- inner_join(LC_meta, biomass_set, by = "ID")

biomass_set <- biomass_set %>% mutate(diam_class = case_when(
  DBH > 0 & DBH <= 20 ~20,
  DBH > 15 & DBH <=40 ~40,
  DBH > 30 & DBH <=60 ~60

  
))

          

biomass_set <- filter(biomass_set, !event_short == "CT1")
str(biomass_set)
biomass_set$DBH <- as.numeric(biomass_set$DBH)
biomass_set$diam_class <- as.factor(biomass_set$diam_class)

biomass_set <- biomass_set %>% group_by(event_short,diam_class)

print(group_data(biomass_set),n=58)






drymass_subset <- subset(biomass_set, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

#how many trees present in each event
drymass_subset <- drymass_subset %>% group_by(event_short)
group_data(drymass_subset)

#13-15E - 31
#16-20 - 31
#2H - 30
#4A - 19
#5A - 31
#5C - 31
#8-9D - 26
#CT3 - 32
#total 

ggplot((drymass_subset), aes(x=Diameter))+
  geom_histogram(bins = 20)

ggplot(subset(drymass_subset, event_short == "13-15E"), aes(x=Diameter))+
         geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "16-20"), aes(x=Diameter))+
  geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "2H"), aes(x=Diameter))+
  geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "4A"), aes(x=Diameter))+
  geom_histogram(bins = 10)

ggplot(subset(drymass_subset, event_short == "5A"), aes(x=Diameter))+
  geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "5C"), aes(x=Diameter))+
  geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "8-9D"), aes(x=Diameter))+
  geom_histogram(bins = 15)

ggplot(subset(drymass_subset, event_short == "CT3"), aes(x=Diameter))+
  geom_histogram(bins = 15)


#Power analysis for volume index


