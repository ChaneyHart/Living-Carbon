#Biomass sampling
library(readxl)
#install.packages("fabricatr")
library(fabricatr)
install.packages("sampling")
library(purrr)
library(tidyr)

growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
#growth_dat$Vol_index <- 0.333*((3.14*((growth_dat$DBH/2)^2))*growth_dat$Height)

growth_dat <- growth_dat %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~"control", 
  event_short == "CT3" ~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

growth_dat <- growth_dat %>% mutate(section = case_when(
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


#Split each event into terciles

biomass_set <- growth_dat %>% group_by(event_short) %>% mutate(
  strata = ntile(V801,3)
)
#Check this worked for each event
biomass_set_5A <- subset(biomass_set, event_short == "5A")
ggplot(biomass_set_5A,aes(x=strata,y=V801))+
  geom_point()
biomass_set_4A <- subset(biomass_set, event_short == "4A")
ggplot(biomass_set_4A,aes(x=strata,y=V801))+
  geom_point()
biomass_set_13_15E <- subset(biomass_set, event_short == "13-15E")
ggplot(biomass_set_13_15E,aes(x=strata,y=V801))+
  geom_point()
biomass_set_16_20 <- subset(biomass_set, event_short == "16-20")
ggplot(biomass_set_16_20,aes(x=strata,y=V801))+
  geom_point()
biomass_set_8_9D <- subset(biomass_set, event_short == "8-9D")
ggplot(biomass_set_8_9D,aes(x=strata,y=V801))+
  geom_point()
biomass_set_CT3 <- subset(biomass_set, event_short == "CT3")
ggplot(biomass_set_CT3,aes(x=strata,y=V801))+
  geom_point()
#see also that the amount of variation is different for each strata

#Select only events of interest
biomass_subset <- subset(biomass_set, event_short == "13-15E" | event_short == "2H" | event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

# add info about the total number of trees in each event (N)
biomass_subset <- biomass_subset %>% mutate(N= case_when(
  event_short == "13-15E"~"32",
  event_short == "16-20"~"33",
  event_short == "5A"~"34",
  event_short == "4A"~"19",
  event_short == "2H"~ "29",
  event_short == "5C"~ "33",
  event_short == "8-9D"~"29",
  event_short == "CT3"~"34"
))



biomass_subset$N <- as.numeric(biomass_subset$N)
#for each event, we want 50% of the tree samples
biomass_subset$n_final <- biomass_subset$N*0.5

#brute force the denominator of Neyman's allocation formula for each event
biomass_subset <- biomass_subset %>% mutate(denom= case_when(
  event_short == "13-15E"~1159.79,
  event_short == "16-20"~1429.84,
  event_short == "5A"~2860.11,
  event_short == "4A"~2170.89,
  event_short == "2H"~ 846.708,
  event_short == "5C"~ 1461.4,
  event_short == "8-9D"~1073.33,
  event_short == "CT3"~724.72
))

#determine size and variance within each strata of each event. Use this to estimate how samples should be allocated to each strata via Neyman allocation
#larger and/or more variable strata get more samples
biomass_summary <- biomass_subset %>% group_by(event_short,strata) %>% summarize(
  ID = ID,
  V = V801,
  n_strata = n(),
  N = N,
  w_strata = n_strata/N,
  sd = sd(V801),
  n_neymans = (n_final)*((w_strata*sd)/denom),
  n_equal = n_final/3,
  n_prop = n_final*w_strata,
  tier = tier,
  block = block,
  construct = construct,
  construct2 = construct2
)


# Take stratified sample.

#make sure data is grouped correctly
print(group_data(biomass_summary),n=24)

nested_neymans <- biomass_summary %>%
  group_by(event_short,strata) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = c(3,3,10,2,3,11,3,4,6,2,2,6,3,3,11,2,3,11,4,2,8,4,4,9))

nested_neymans$event_short <- as.factor(nested_neymans$event_short)
nested_neymans$strata <- as.factor(nested_neymans$strata)
print(nested_neymans,n=24)

sampled_neymans <- nested_neymans %>%
    mutate(samp = map2(data, n, sample_n))

sampled_neymans_list <- sampled_neymans %>%
  select(-data) %>%
  unnest(samp)


sampled_neymans_list <- subset(sampled_neymans_list, event_short != "2H" & event_short != "5C")

sample_list_priority <- subset(sampled_neymans_list, event_short == "5A" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D")

sample_list_2nd_priority <- subset(sampled_neymans_list, event_short == "13-15E" | event_short == "CT3")


