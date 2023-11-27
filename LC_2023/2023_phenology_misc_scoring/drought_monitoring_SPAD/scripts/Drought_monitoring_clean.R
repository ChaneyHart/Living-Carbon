#LC growth analysis update - Sep 21 2023


library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)

#read in data

drought_scores_SPAD <- read.csv("LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/Drought_Scores_SPAD_compiled.csv")
growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
drought_scores_SPAD <- inner_join(drought_scores_SPAD,growth)

#seperate so NAs in one do not affect the other
drought_scores_2023 <- subset(drought_scores_SPAD, select = c("row","column","event_short","construct","construct2","ID","DS_1","DS_18","DS_32","DS_55"))
SPAD_2023 <- subset(drought_scores_SPAD, select = c("row","column","event_short","construct","construct2","ID","SPAD_18","SPAD_32","SPAD_55"))

#analyzing drought scores
drought_scores_2023 <- na.omit(drought_scores_2023)

#drought score is subjective but we will assume it does not go backwards

drought_scores_2023$DS_18 <- as.numeric(drought_scores_2023$DS_18)
drought_scores_2023$DS_32 <- as.numeric(drought_scores_2023$DS_32)
drought_scores_2023$DS_55 <- as.numeric(drought_scores_2023$DS_55)

drought_scores_2023$delta1 <- drought_scores_2023$DS_18 - drought_scores_2023$DS_1

drought_delta_1 <- subset(drought_scores_2023, delta1 < 0)
delta1 <- c(row.names(drought_delta_1))

for (i in delta1) {
  print(i)
  drought_scores_2023[delta1, 8] <- drought_scores_2023[delta1, 7]
}


drought_scores_2023$delta2 <- drought_scores_2023$DS_32 - drought_scores_2023$DS_18

drought_delta_2 <- subset(drought_scores_2023, delta2 < 0)
delta2 <- c(row.names(drought_delta_2))

for (i in delta2) {
  print(i)
  drought_scores_2023[delta2, 9] <- drought_scores_2023[delta2, 8]
}  



drought_scores_2023$delta3 <- drought_scores_2023$DS_55 - drought_scores_2023$DS_32

drought_delta_3 <- subset(drought_scores_2023, delta3 < 0)
delta3 <- c(row.names(drought_delta_3))

for (i in delta3) {
  print(i)
  drought_scores_2023[delta3, 10] <- drought_scores_2023[delta3, 9]
}  

write.csv(file = "LC_2023/2023_phenology_misc_scoring/Drought_scores_clean_2023.csv",drought_scores_2023)


###########################
### summarize proportions of scores for each event


######CT3

drought_CT3 <- subset(drought_scores_2023, event_short == "CT3")

drought_CT3_summary <- drought_CT3 %>% summarize(
  day01_0 = count(subset(drought_CT3, DS_1 == "0")),
  day01_1 = count(subset(drought_CT3, DS_1 == "1")),
  day01_2 = count(subset(drought_CT3, DS_1 == "2")),
  day01_3 = count(subset(drought_CT3, DS_1 == "3")),
  day01_4 = count(subset(drought_CT3, DS_1 == "4")),
  day01_5 = count(subset(drought_CT3, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_CT3, DS_18 == "0")),
  day18_1 = count(subset(drought_CT3, DS_18 == "1")),
  day18_2 = count(subset(drought_CT3, DS_18 == "2")),
  day18_3 = count(subset(drought_CT3, DS_18 == "3")),
  day18_4 = count(subset(drought_CT3, DS_18 == "4")),
  day18_5 = count(subset(drought_CT3, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_CT3, DS_32 == "0")),
  day32_1 = count(subset(drought_CT3, DS_32 == "1")),
  day32_2 = count(subset(drought_CT3, DS_32 == "2")),
  day32_3 = count(subset(drought_CT3, DS_32 == "3")),
  day32_4 = count(subset(drought_CT3, DS_32 == "4")),
  day32_5 = count(subset(drought_CT3, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_CT3, DS_55 == "0")),
  day55_1 = count(subset(drought_CT3, DS_55 == "1")),
  day55_2 = count(subset(drought_CT3, DS_55 == "2")),
  day55_3 = count(subset(drought_CT3, DS_55 == "3")),
  day55_4 = count(subset(drought_CT3, DS_55 == "4")),
  day55_5 = count(subset(drought_CT3, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_CT3_summary2 <- drought_CT3_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_CT3_summary2_0 <- pivot_longer(drought_CT3_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_CT3_summary2_0 <- drought_CT3_summary2_0[,21:22]
drought_CT3_summary2_0$event <- "CT3"
drought_CT3_summary2_0$construct <- "Control"
drought_CT3_summary2_0$Day <- substr(drought_CT3_summary2_0$Day,4,5)
drought_CT3_summary2_0 <- as.data.frame(as.matrix(drought_CT3_summary2_0[,])) 

#1
drought_CT3_summary2_1 <- pivot_longer(drought_CT3_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_CT3_summary2_1 <- drought_CT3_summary2_1[,21:22]
drought_CT3_summary2_1$event <- "CT3"
drought_CT3_summary2_1$construct <- "Control"
drought_CT3_summary2_1$Day <- substr(drought_CT3_summary2_1$Day,4,5)
drought_CT3_summary2_1 <- as.data.frame(as.matrix(drought_CT3_summary2_1[,])) 

#2
drought_CT3_summary2_2 <- pivot_longer(drought_CT3_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_CT3_summary2_2 <- drought_CT3_summary2_2[,21:22]
drought_CT3_summary2_2$event <- "CT3"
drought_CT3_summary2_2$construct <- "Control"
drought_CT3_summary2_2$Day <- substr(drought_CT3_summary2_2$Day,4,5)
drought_CT3_summary2_2 <- as.data.frame(as.matrix(drought_CT3_summary2_2[,])) 

#3
drought_CT3_summary2_3 <- pivot_longer(drought_CT3_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_CT3_summary2_3 <- drought_CT3_summary2_3[,21:22]
drought_CT3_summary2_3$event <- "CT3"
drought_CT3_summary2_3$construct <- "Control"
drought_CT3_summary2_3$Day <- substr(drought_CT3_summary2_3$Day,4,5)
drought_CT3_summary2_3 <- as.data.frame(as.matrix(drought_CT3_summary2_3[,])) 
#4
drought_CT3_summary2_4 <- pivot_longer(drought_CT3_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_CT3_summary2_4 <- drought_CT3_summary2_4[,21:22]
drought_CT3_summary2_4$event <- "CT3"
drought_CT3_summary2_4$construct <- "Control"
drought_CT3_summary2_4$Day <- substr(drought_CT3_summary2_4$Day,4,5)
drought_CT3_summary2_4 <- as.data.frame(as.matrix(drought_CT3_summary2_4[,])) 

#5
drought_CT3_summary2_5 <- pivot_longer(drought_CT3_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_CT3_summary2_5 <- drought_CT3_summary2_5[,21:22]
drought_CT3_summary2_5$event <- "CT3"
drought_CT3_summary2_5$construct <- "Control"
drought_CT3_summary2_5$Day <- substr(drought_CT3_summary2_5$Day,4,5)
drought_CT3_summary2_5 <- as.data.frame(as.matrix(drought_CT3_summary2_5[,])) 

drought_CT3_summary_3 <- inner_join(drought_CT3_summary2_0, drought_CT3_summary2_1, by = c("Day","event","construct"))
drought_CT3_summary_3 <- inner_join(drought_CT3_summary_3, drought_CT3_summary2_2, by = c("Day","event","construct"))
drought_CT3_summary_3 <- inner_join(drought_CT3_summary_3, drought_CT3_summary2_3, by = c("Day","event","construct"))
drought_CT3_summary_3 <- inner_join(drought_CT3_summary_3, drought_CT3_summary2_4, by = c("Day","event","construct"))
drought_CT3_summary_3 <- inner_join(drought_CT3_summary_3, drought_CT3_summary2_5, by = c("Day","event","construct"))

#####
# 16-20

drought_16_20 <- subset(drought_scores_2023, event_short == "16-20")

drought_16_20_summary <- drought_16_20 %>% summarize(
  day01_0 = count(subset(drought_16_20, DS_1 == "0")),
  day01_1 = count(subset(drought_16_20, DS_1 == "1")),
  day01_2 = count(subset(drought_16_20, DS_1 == "2")),
  day01_3 = count(subset(drought_16_20, DS_1 == "3")),
  day01_4 = count(subset(drought_16_20, DS_1 == "4")),
  day01_5 = count(subset(drought_16_20, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_16_20, DS_18 == "0")),
  day18_1 = count(subset(drought_16_20, DS_18 == "1")),
  day18_2 = count(subset(drought_16_20, DS_18 == "2")),
  day18_3 = count(subset(drought_16_20, DS_18 == "3")),
  day18_4 = count(subset(drought_16_20, DS_18 == "4")),
  day18_5 = count(subset(drought_16_20, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_16_20, DS_32 == "0")),
  day32_1 = count(subset(drought_16_20, DS_32 == "1")),
  day32_2 = count(subset(drought_16_20, DS_32 == "2")),
  day32_3 = count(subset(drought_16_20, DS_32 == "3")),
  day32_4 = count(subset(drought_16_20, DS_32 == "4")),
  day32_5 = count(subset(drought_16_20, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_16_20, DS_55 == "0")),
  day55_1 = count(subset(drought_16_20, DS_55 == "1")),
  day55_2 = count(subset(drought_16_20, DS_55 == "2")),
  day55_3 = count(subset(drought_16_20, DS_55 == "3")),
  day55_4 = count(subset(drought_16_20, DS_55 == "4")),
  day55_5 = count(subset(drought_16_20, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_16_20_summary2 <- drought_16_20_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_16_20_summary2_0 <- pivot_longer(drought_16_20_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_16_20_summary2_0 <- drought_16_20_summary2_0[,21:22]
drought_16_20_summary2_0$event <- "16-20"
drought_16_20_summary2_0$construct <- "Escape"
drought_16_20_summary2_0$Day <- substr(drought_16_20_summary2_0$Day,4,5)
drought_16_20_summary2_0 <- as.data.frame(as.matrix(drought_16_20_summary2_0[,])) 

#1
drought_16_20_summary2_1 <- pivot_longer(drought_16_20_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_16_20_summary2_1 <- drought_16_20_summary2_1[,21:22]
drought_16_20_summary2_1$event <- "16-20"
drought_16_20_summary2_1$construct <- "Escape"
drought_16_20_summary2_1$Day <- substr(drought_16_20_summary2_1$Day,4,5)
drought_16_20_summary2_1 <- as.data.frame(as.matrix(drought_16_20_summary2_1[,])) 

#2
drought_16_20_summary2_2 <- pivot_longer(drought_16_20_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_16_20_summary2_2 <- drought_16_20_summary2_2[,21:22]
drought_16_20_summary2_2$event <- "16-20"
drought_16_20_summary2_2$construct <- "Escape"
drought_16_20_summary2_2$Day <- substr(drought_16_20_summary2_2$Day,4,5)
drought_16_20_summary2_2 <- as.data.frame(as.matrix(drought_16_20_summary2_2[,])) 

#3
drought_16_20_summary2_3 <- pivot_longer(drought_16_20_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_16_20_summary2_3 <- drought_16_20_summary2_3[,21:22]
drought_16_20_summary2_3$event <- "16-20"
drought_16_20_summary2_3$construct <- "Escape"
drought_16_20_summary2_3$Day <- substr(drought_16_20_summary2_3$Day,4,5)
drought_16_20_summary2_3 <- as.data.frame(as.matrix(drought_16_20_summary2_3[,])) 
#4
drought_16_20_summary2_4 <- pivot_longer(drought_16_20_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_16_20_summary2_4 <- drought_16_20_summary2_4[,21:22]
drought_16_20_summary2_4$event <- "16-20"
drought_16_20_summary2_4$construct <- "Escape"
drought_16_20_summary2_4$Day <- substr(drought_16_20_summary2_4$Day,4,5)
drought_16_20_summary2_4 <- as.data.frame(as.matrix(drought_16_20_summary2_4[,])) 

#5
drought_16_20_summary2_5 <- pivot_longer(drought_16_20_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_16_20_summary2_5 <- drought_16_20_summary2_5[,21:22]
drought_16_20_summary2_5$event <- "16-20"
drought_16_20_summary2_5$construct <- "Escape"
drought_16_20_summary2_5$Day <- substr(drought_16_20_summary2_5$Day,4,5)
drought_16_20_summary2_5 <- as.data.frame(as.matrix(drought_16_20_summary2_5[,])) 

drought_16_20_summary_3 <- inner_join(drought_16_20_summary2_0, drought_16_20_summary2_1, by = c("Day","event","construct"))
drought_16_20_summary_3 <- inner_join(drought_16_20_summary_3, drought_16_20_summary2_2, by = c("Day","event","construct"))
drought_16_20_summary_3 <- inner_join(drought_16_20_summary_3, drought_16_20_summary2_3, by = c("Day","event","construct"))
drought_16_20_summary_3 <- inner_join(drought_16_20_summary_3, drought_16_20_summary2_4, by = c("Day","event","construct"))
drought_16_20_summary_3 <- inner_join(drought_16_20_summary_3, drought_16_20_summary2_5, by = c("Day","event","construct"))

####
# 8-9D

drought_8_9D <- subset(drought_scores_2023, event_short == "8-9D")

drought_8_9D_summary <- drought_8_9D %>% summarize(
  day01_0 = count(subset(drought_8_9D, DS_1 == "0")),
  day01_1 = count(subset(drought_8_9D, DS_1 == "1")),
  day01_2 = count(subset(drought_8_9D, DS_1 == "2")),
  day01_3 = count(subset(drought_8_9D, DS_1 == "3")),
  day01_4 = count(subset(drought_8_9D, DS_1 == "4")),
  day01_5 = count(subset(drought_8_9D, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_8_9D, DS_18 == "0")),
  day18_1 = count(subset(drought_8_9D, DS_18 == "1")),
  day18_2 = count(subset(drought_8_9D, DS_18 == "2")),
  day18_3 = count(subset(drought_8_9D, DS_18 == "3")),
  day18_4 = count(subset(drought_8_9D, DS_18 == "4")),
  day18_5 = count(subset(drought_8_9D, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_8_9D, DS_32 == "0")),
  day32_1 = count(subset(drought_8_9D, DS_32 == "1")),
  day32_2 = count(subset(drought_8_9D, DS_32 == "2")),
  day32_3 = count(subset(drought_8_9D, DS_32 == "3")),
  day32_4 = count(subset(drought_8_9D, DS_32 == "4")),
  day32_5 = count(subset(drought_8_9D, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_8_9D, DS_55 == "0")),
  day55_1 = count(subset(drought_8_9D, DS_55 == "1")),
  day55_2 = count(subset(drought_8_9D, DS_55 == "2")),
  day55_3 = count(subset(drought_8_9D, DS_55 == "3")),
  day55_4 = count(subset(drought_8_9D, DS_55 == "4")),
  day55_5 = count(subset(drought_8_9D, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_8_9D_summary2 <- drought_8_9D_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_8_9D_summary2_0 <- pivot_longer(drought_8_9D_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_8_9D_summary2_0 <- drought_8_9D_summary2_0[,21:22]
drought_8_9D_summary2_0$event <- "8-9D"
drought_8_9D_summary2_0$construct <- "Escape"
drought_8_9D_summary2_0$Day <- substr(drought_8_9D_summary2_0$Day,4,5)
drought_8_9D_summary2_0 <- as.data.frame(as.matrix(drought_8_9D_summary2_0[,])) 

#1
drought_8_9D_summary2_1 <- pivot_longer(drought_8_9D_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_8_9D_summary2_1 <- drought_8_9D_summary2_1[,21:22]
drought_8_9D_summary2_1$event <- "8-9D"
drought_8_9D_summary2_1$construct <- "Escape"
drought_8_9D_summary2_1$Day <- substr(drought_8_9D_summary2_1$Day,4,5)
drought_8_9D_summary2_1 <- as.data.frame(as.matrix(drought_8_9D_summary2_1[,])) 

#2
drought_8_9D_summary2_2 <- pivot_longer(drought_8_9D_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_8_9D_summary2_2 <- drought_8_9D_summary2_2[,21:22]
drought_8_9D_summary2_2$event <- "8-9D"
drought_8_9D_summary2_2$construct <- "Escape"
drought_8_9D_summary2_2$Day <- substr(drought_8_9D_summary2_2$Day,4,5)
drought_8_9D_summary2_2 <- as.data.frame(as.matrix(drought_8_9D_summary2_2[,])) 

#3
drought_8_9D_summary2_3 <- pivot_longer(drought_8_9D_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_8_9D_summary2_3 <- drought_8_9D_summary2_3[,21:22]
drought_8_9D_summary2_3$event <- "8-9D"
drought_8_9D_summary2_3$construct <- "Escape"
drought_8_9D_summary2_3$Day <- substr(drought_8_9D_summary2_3$Day,4,5)
drought_8_9D_summary2_3 <- as.data.frame(as.matrix(drought_8_9D_summary2_3[,])) 
#4
drought_8_9D_summary2_4 <- pivot_longer(drought_8_9D_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_8_9D_summary2_4 <- drought_8_9D_summary2_4[,21:22]
drought_8_9D_summary2_4$event <- "8-9D"
drought_8_9D_summary2_4$construct <- "Escape"
drought_8_9D_summary2_4$Day <- substr(drought_8_9D_summary2_4$Day,4,5)
drought_8_9D_summary2_4 <- as.data.frame(as.matrix(drought_8_9D_summary2_4[,])) 

#5
drought_8_9D_summary2_5 <- pivot_longer(drought_8_9D_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_8_9D_summary2_5 <- drought_8_9D_summary2_5[,21:22]
drought_8_9D_summary2_5$event <- "8-9D"
drought_8_9D_summary2_5$construct <- "Escape"
drought_8_9D_summary2_5$Day <- substr(drought_8_9D_summary2_5$Day,4,5)
drought_8_9D_summary2_5 <- as.data.frame(as.matrix(drought_8_9D_summary2_5[,])) 

drought_8_9D_summary_3 <- inner_join(drought_8_9D_summary2_0, drought_8_9D_summary2_1, by = c("Day","event","construct"))
drought_8_9D_summary_3 <- inner_join(drought_8_9D_summary_3, drought_8_9D_summary2_2, by = c("Day","event","construct"))
drought_8_9D_summary_3 <- inner_join(drought_8_9D_summary_3, drought_8_9D_summary2_3, by = c("Day","event","construct"))
drought_8_9D_summary_3 <- inner_join(drought_8_9D_summary_3, drought_8_9D_summary2_4, by = c("Day","event","construct"))
drought_8_9D_summary_3 <- inner_join(drought_8_9D_summary_3, drought_8_9D_summary2_5, by = c("Day","event","construct"))

####
# 5A

drought_5A <- subset(drought_scores_2023, event_short == "5A")

drought_5A_summary <- drought_5A %>% summarize(
  day01_0 = count(subset(drought_5A, DS_1 == "0")),
  day01_1 = count(subset(drought_5A, DS_1 == "1")),
  day01_2 = count(subset(drought_5A, DS_1 == "2")),
  day01_3 = count(subset(drought_5A, DS_1 == "3")),
  day01_4 = count(subset(drought_5A, DS_1 == "4")),
  day01_5 = count(subset(drought_5A, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_5A, DS_18 == "0")),
  day18_1 = count(subset(drought_5A, DS_18 == "1")),
  day18_2 = count(subset(drought_5A, DS_18 == "2")),
  day18_3 = count(subset(drought_5A, DS_18 == "3")),
  day18_4 = count(subset(drought_5A, DS_18 == "4")),
  day18_5 = count(subset(drought_5A, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_5A, DS_32 == "0")),
  day32_1 = count(subset(drought_5A, DS_32 == "1")),
  day32_2 = count(subset(drought_5A, DS_32 == "2")),
  day32_3 = count(subset(drought_5A, DS_32 == "3")),
  day32_4 = count(subset(drought_5A, DS_32 == "4")),
  day32_5 = count(subset(drought_5A, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_5A, DS_55 == "0")),
  day55_1 = count(subset(drought_5A, DS_55 == "1")),
  day55_2 = count(subset(drought_5A, DS_55 == "2")),
  day55_3 = count(subset(drought_5A, DS_55 == "3")),
  day55_4 = count(subset(drought_5A, DS_55 == "4")),
  day55_5 = count(subset(drought_5A, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_5A_summary2 <- drought_5A_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_5A_summary2_0 <- pivot_longer(drought_5A_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_5A_summary2_0 <- drought_5A_summary2_0[,21:22]
drought_5A_summary2_0$event <- "5A"
drought_5A_summary2_0$construct <- "Transgenic"
drought_5A_summary2_0$Day <- substr(drought_5A_summary2_0$Day,4,5)
drought_5A_summary2_0 <- as.data.frame(as.matrix(drought_5A_summary2_0[,])) 

#1
drought_5A_summary2_1 <- pivot_longer(drought_5A_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_5A_summary2_1 <- drought_5A_summary2_1[,21:22]
drought_5A_summary2_1$event <- "5A"
drought_5A_summary2_1$construct <- "Transgenic"
drought_5A_summary2_1$Day <- substr(drought_5A_summary2_1$Day,4,5)
drought_5A_summary2_1 <- as.data.frame(as.matrix(drought_5A_summary2_1[,])) 

#2
drought_5A_summary2_2 <- pivot_longer(drought_5A_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_5A_summary2_2 <- drought_5A_summary2_2[,21:22]
drought_5A_summary2_2$event <- "5A"
drought_5A_summary2_2$construct <- "Transgenic"
drought_5A_summary2_2$Day <- substr(drought_5A_summary2_2$Day,4,5)
drought_5A_summary2_2 <- as.data.frame(as.matrix(drought_5A_summary2_2[,])) 

#3
drought_5A_summary2_3 <- pivot_longer(drought_5A_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_5A_summary2_3 <- drought_5A_summary2_3[,21:22]
drought_5A_summary2_3$event <- "5A"
drought_5A_summary2_3$construct <- "Transgenic"
drought_5A_summary2_3$Day <- substr(drought_5A_summary2_3$Day,4,5)
drought_5A_summary2_3 <- as.data.frame(as.matrix(drought_5A_summary2_3[,])) 
#4
drought_5A_summary2_4 <- pivot_longer(drought_5A_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_5A_summary2_4 <- drought_5A_summary2_4[,21:22]
drought_5A_summary2_4$event <- "5A"
drought_5A_summary2_4$construct <- "Transgenic"
drought_5A_summary2_4$Day <- substr(drought_5A_summary2_4$Day,4,5)
drought_5A_summary2_4 <- as.data.frame(as.matrix(drought_5A_summary2_4[,])) 

#5
drought_5A_summary2_5 <- pivot_longer(drought_5A_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_5A_summary2_5 <- drought_5A_summary2_5[,21:22]
drought_5A_summary2_5$event <- "5A"
drought_5A_summary2_5$construct <- "Transgenic"
drought_5A_summary2_5$Day <- substr(drought_5A_summary2_5$Day,4,5)
drought_5A_summary2_5 <- as.data.frame(as.matrix(drought_5A_summary2_5[,])) 

drought_5A_summary_3 <- inner_join(drought_5A_summary2_0, drought_5A_summary2_1, by = c("Day","event","construct"))
drought_5A_summary_3 <- inner_join(drought_5A_summary_3, drought_5A_summary2_2, by = c("Day","event","construct"))
drought_5A_summary_3 <- inner_join(drought_5A_summary_3, drought_5A_summary2_3, by = c("Day","event","construct"))
drought_5A_summary_3 <- inner_join(drought_5A_summary_3, drought_5A_summary2_4, by = c("Day","event","construct"))
drought_5A_summary_3 <- inner_join(drought_5A_summary_3, drought_5A_summary2_5, by = c("Day","event","construct"))

####
#4A

drought_4A <- subset(drought_scores_2023, event_short == "4A")

drought_4A_summary <- drought_4A %>% summarize(
  day01_0 = count(subset(drought_4A, DS_1 == "0")),
  day01_1 = count(subset(drought_4A, DS_1 == "1")),
  day01_2 = count(subset(drought_4A, DS_1 == "2")),
  day01_3 = count(subset(drought_4A, DS_1 == "3")),
  day01_4 = count(subset(drought_4A, DS_1 == "4")),
  day01_5 = count(subset(drought_4A, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_4A, DS_18 == "0")),
  day18_1 = count(subset(drought_4A, DS_18 == "1")),
  day18_2 = count(subset(drought_4A, DS_18 == "2")),
  day18_3 = count(subset(drought_4A, DS_18 == "3")),
  day18_4 = count(subset(drought_4A, DS_18 == "4")),
  day18_5 = count(subset(drought_4A, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_4A, DS_32 == "0")),
  day32_1 = count(subset(drought_4A, DS_32 == "1")),
  day32_2 = count(subset(drought_4A, DS_32 == "2")),
  day32_3 = count(subset(drought_4A, DS_32 == "3")),
  day32_4 = count(subset(drought_4A, DS_32 == "4")),
  day32_5 = count(subset(drought_4A, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_4A, DS_55 == "0")),
  day55_1 = count(subset(drought_4A, DS_55 == "1")),
  day55_2 = count(subset(drought_4A, DS_55 == "2")),
  day55_3 = count(subset(drought_4A, DS_55 == "3")),
  day55_4 = count(subset(drought_4A, DS_55 == "4")),
  day55_5 = count(subset(drought_4A, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_4A_summary2 <- drought_4A_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_4A_summary2_0 <- pivot_longer(drought_4A_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_4A_summary2_0 <- drought_4A_summary2_0[,21:22]
drought_4A_summary2_0$event <- "4A"
drought_4A_summary2_0$construct <- "Transgenic"
drought_4A_summary2_0$Day <- substr(drought_4A_summary2_0$Day,4,5)
drought_4A_summary2_0 <- as.data.frame(as.matrix(drought_4A_summary2_0[,])) 

#1
drought_4A_summary2_1 <- pivot_longer(drought_4A_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_4A_summary2_1 <- drought_4A_summary2_1[,21:22]
drought_4A_summary2_1$event <- "4A"
drought_4A_summary2_1$construct <- "Transgenic"
drought_4A_summary2_1$Day <- substr(drought_4A_summary2_1$Day,4,5)
drought_4A_summary2_1 <- as.data.frame(as.matrix(drought_4A_summary2_1[,])) 

#2
drought_4A_summary2_2 <- pivot_longer(drought_4A_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_4A_summary2_2 <- drought_4A_summary2_2[,21:22]
drought_4A_summary2_2$event <- "4A"
drought_4A_summary2_2$construct <- "Transgenic"
drought_4A_summary2_2$Day <- substr(drought_4A_summary2_2$Day,4,5)
drought_4A_summary2_2 <- as.data.frame(as.matrix(drought_4A_summary2_2[,])) 

#3
drought_4A_summary2_3 <- pivot_longer(drought_4A_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_4A_summary2_3 <- drought_4A_summary2_3[,21:22]
drought_4A_summary2_3$event <- "4A"
drought_4A_summary2_3$construct <- "Transgenic"
drought_4A_summary2_3$Day <- substr(drought_4A_summary2_3$Day,4,5)
drought_4A_summary2_3 <- as.data.frame(as.matrix(drought_4A_summary2_3[,])) 
#4
drought_4A_summary2_4 <- pivot_longer(drought_4A_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_4A_summary2_4 <- drought_4A_summary2_4[,21:22]
drought_4A_summary2_4$event <- "4A"
drought_4A_summary2_4$construct <- "Transgenic"
drought_4A_summary2_4$Day <- substr(drought_4A_summary2_4$Day,4,5)
drought_4A_summary2_4 <- as.data.frame(as.matrix(drought_4A_summary2_4[,])) 

#5
drought_4A_summary2_5 <- pivot_longer(drought_4A_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_4A_summary2_5 <- drought_4A_summary2_5[,21:22]
drought_4A_summary2_5$event <- "4A"
drought_4A_summary2_5$construct <- "Transgenic"
drought_4A_summary2_5$Day <- substr(drought_4A_summary2_5$Day,4,5)
drought_4A_summary2_5 <- as.data.frame(as.matrix(drought_4A_summary2_5[,])) 

drought_4A_summary_3 <- inner_join(drought_4A_summary2_0, drought_4A_summary2_1, by = c("Day","event","construct"))
drought_4A_summary_3 <- inner_join(drought_4A_summary_3, drought_4A_summary2_2, by = c("Day","event","construct"))
drought_4A_summary_3 <- inner_join(drought_4A_summary_3, drought_4A_summary2_3, by = c("Day","event","construct"))
drought_4A_summary_3 <- inner_join(drought_4A_summary_3, drought_4A_summary2_4, by = c("Day","event","construct"))
drought_4A_summary_3 <- inner_join(drought_4A_summary_3, drought_4A_summary2_5, by = c("Day","event","construct"))

####
#13-15E

drought_13_15E <- subset(drought_scores_2023, event_short == "13-15E")

drought_13_15E_summary <- drought_13_15E %>% summarize(
  day01_0 = count(subset(drought_13_15E, DS_1 == "0")),
  day01_1 = count(subset(drought_13_15E, DS_1 == "1")),
  day01_2 = count(subset(drought_13_15E, DS_1 == "2")),
  day01_3 = count(subset(drought_13_15E, DS_1 == "3")),
  day01_4 = count(subset(drought_13_15E, DS_1 == "4")),
  day01_5 = count(subset(drought_13_15E, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_13_15E, DS_18 == "0")),
  day18_1 = count(subset(drought_13_15E, DS_18 == "1")),
  day18_2 = count(subset(drought_13_15E, DS_18 == "2")),
  day18_3 = count(subset(drought_13_15E, DS_18 == "3")),
  day18_4 = count(subset(drought_13_15E, DS_18 == "4")),
  day18_5 = count(subset(drought_13_15E, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_13_15E, DS_32 == "0")),
  day32_1 = count(subset(drought_13_15E, DS_32 == "1")),
  day32_2 = count(subset(drought_13_15E, DS_32 == "2")),
  day32_3 = count(subset(drought_13_15E, DS_32 == "3")),
  day32_4 = count(subset(drought_13_15E, DS_32 == "4")),
  day32_5 = count(subset(drought_13_15E, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_13_15E, DS_55 == "0")),
  day55_1 = count(subset(drought_13_15E, DS_55 == "1")),
  day55_2 = count(subset(drought_13_15E, DS_55 == "2")),
  day55_3 = count(subset(drought_13_15E, DS_55 == "3")),
  day55_4 = count(subset(drought_13_15E, DS_55 == "4")),
  day55_5 = count(subset(drought_13_15E, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_13_15E_summary2 <- drought_13_15E_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_13_15E_summary2_0 <- pivot_longer(drought_13_15E_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_13_15E_summary2_0 <- drought_13_15E_summary2_0[,21:22]
drought_13_15E_summary2_0$event <- "13-15E"
drought_13_15E_summary2_0$construct <- "Transgenic"
drought_13_15E_summary2_0$Day <- substr(drought_13_15E_summary2_0$Day,4,5)
drought_13_15E_summary2_0 <- as.data.frame(as.matrix(drought_13_15E_summary2_0[,])) 

#1
drought_13_15E_summary2_1 <- pivot_longer(drought_13_15E_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_13_15E_summary2_1 <- drought_13_15E_summary2_1[,21:22]
drought_13_15E_summary2_1$event <- "13-15E"
drought_13_15E_summary2_1$construct <- "Transgenic"
drought_13_15E_summary2_1$Day <- substr(drought_13_15E_summary2_1$Day,4,5)
drought_13_15E_summary2_1 <- as.data.frame(as.matrix(drought_13_15E_summary2_1[,])) 

#2
drought_13_15E_summary2_2 <- pivot_longer(drought_13_15E_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_13_15E_summary2_2 <- drought_13_15E_summary2_2[,21:22]
drought_13_15E_summary2_2$event <- "13-15E"
drought_13_15E_summary2_2$construct <- "Transgenic"
drought_13_15E_summary2_2$Day <- substr(drought_13_15E_summary2_2$Day,4,5)
drought_13_15E_summary2_2 <- as.data.frame(as.matrix(drought_13_15E_summary2_2[,])) 

#3
drought_13_15E_summary2_3 <- pivot_longer(drought_13_15E_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_13_15E_summary2_3 <- drought_13_15E_summary2_3[,21:22]
drought_13_15E_summary2_3$event <- "13-15E"
drought_13_15E_summary2_3$construct <- "Transgenic"
drought_13_15E_summary2_3$Day <- substr(drought_13_15E_summary2_3$Day,4,5)
drought_13_15E_summary2_3 <- as.data.frame(as.matrix(drought_13_15E_summary2_3[,])) 
#4
drought_13_15E_summary2_4 <- pivot_longer(drought_13_15E_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_13_15E_summary2_4 <- drought_13_15E_summary2_4[,21:22]
drought_13_15E_summary2_4$event <- "13-15E"
drought_13_15E_summary2_4$construct <- "Transgenic"
drought_13_15E_summary2_4$Day <- substr(drought_13_15E_summary2_4$Day,4,5)
drought_13_15E_summary2_4 <- as.data.frame(as.matrix(drought_13_15E_summary2_4[,])) 

#5
drought_13_15E_summary2_5 <- pivot_longer(drought_13_15E_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_13_15E_summary2_5 <- drought_13_15E_summary2_5[,21:22]
drought_13_15E_summary2_5$event <- "13-15E"
drought_13_15E_summary2_5$construct <- "Transgenic"
drought_13_15E_summary2_5$Day <- substr(drought_13_15E_summary2_5$Day,4,5)
drought_13_15E_summary2_5 <- as.data.frame(as.matrix(drought_13_15E_summary2_5[,])) 

drought_13_15E_summary_3 <- inner_join(drought_13_15E_summary2_0, drought_13_15E_summary2_1, by = c("Day","event","construct"))
drought_13_15E_summary_3 <- inner_join(drought_13_15E_summary_3, drought_13_15E_summary2_2, by = c("Day","event","construct"))
drought_13_15E_summary_3 <- inner_join(drought_13_15E_summary_3, drought_13_15E_summary2_3, by = c("Day","event","construct"))
drought_13_15E_summary_3 <- inner_join(drought_13_15E_summary_3, drought_13_15E_summary2_4, by = c("Day","event","construct"))
drought_13_15E_summary_3 <- inner_join(drought_13_15E_summary_3, drought_13_15E_summary2_5, by = c("Day","event","construct"))

####
#2H

drought_2H <- subset(drought_scores_2023, event_short == "2H")

drought_2H_summary <- drought_2H %>% summarize(
  day01_0 = count(subset(drought_2H, DS_1 == "0")),
  day01_1 = count(subset(drought_2H, DS_1 == "1")),
  day01_2 = count(subset(drought_2H, DS_1 == "2")),
  day01_3 = count(subset(drought_2H, DS_1 == "3")),
  day01_4 = count(subset(drought_2H, DS_1 == "4")),
  day01_5 = count(subset(drought_2H, DS_1 == "5")),
  day01_prop_0 = day01_0/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_1 = day01_1/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_2 = day01_2/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_3 = day01_3/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_4 = day01_4/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day01_prop_5 = day01_5/(sum(day01_0,day01_1,day01_2,day01_3,day01_4,day01_5)),
  day18_0 = count(subset(drought_2H, DS_18 == "0")),
  day18_1 = count(subset(drought_2H, DS_18 == "1")),
  day18_2 = count(subset(drought_2H, DS_18 == "2")),
  day18_3 = count(subset(drought_2H, DS_18 == "3")),
  day18_4 = count(subset(drought_2H, DS_18 == "4")),
  day18_5 = count(subset(drought_2H, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_1 = day18_1/(sum(day18_1,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_2 = day18_2/(sum(day18_2,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_3 = day18_3/(sum(day18_3,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_4 = day18_4/(sum(day18_4,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day18_prop_5 = day18_4/(sum(day18_5,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_2H, DS_32 == "0")),
  day32_1 = count(subset(drought_2H, DS_32 == "1")),
  day32_2 = count(subset(drought_2H, DS_32 == "2")),
  day32_3 = count(subset(drought_2H, DS_32 == "3")),
  day32_4 = count(subset(drought_2H, DS_32 == "4")),
  day32_5 = count(subset(drought_2H, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_1 = day32_1/(sum(day32_1,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_2 = day32_2/(sum(day32_2,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_3 = day32_3/(sum(day32_3,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_4 = day32_4/(sum(day32_4,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day32_prop_5 = day32_5/(sum(day32_5,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_2H, DS_55 == "0")),
  day55_1 = count(subset(drought_2H, DS_55 == "1")),
  day55_2 = count(subset(drought_2H, DS_55 == "2")),
  day55_3 = count(subset(drought_2H, DS_55 == "3")),
  day55_4 = count(subset(drought_2H, DS_55 == "4")),
  day55_5 = count(subset(drought_2H, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_1 = day55_1/(sum(day55_1,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_2 = day55_2/(sum(day55_2,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_3 = day55_3/(sum(day55_3,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_4 = day55_4/(sum(day55_4,day55_1,day55_2,day55_3,day55_4,day55_5)),
  day55_prop_5 = day55_5/(sum(day55_5,day55_1,day55_2,day55_3,day55_4,day55_5)),
)
drought_2H_summary2 <- drought_2H_summary[1,c(7:12,19:24,31:36,43:48)]
#reshape data frame for each proportion
#0
drought_2H_summary2_0 <- pivot_longer(drought_2H_summary2, cols = c(1,7,13,19), names_to = "Day",values_to = "percent_0")
drought_2H_summary2_0 <- drought_2H_summary2_0[,21:22]
drought_2H_summary2_0$event <- "2H"
drought_2H_summary2_0$construct <- "Transgenic"
drought_2H_summary2_0$Day <- substr(drought_2H_summary2_0$Day,4,5)
drought_2H_summary2_0 <- as.data.frame(as.matrix(drought_2H_summary2_0[,])) 

#1
drought_2H_summary2_1 <- pivot_longer(drought_2H_summary2, cols = c(2,8,14,20), names_to = "Day",values_to = "percent_1")
drought_2H_summary2_1 <- drought_2H_summary2_1[,21:22]
drought_2H_summary2_1$event <- "2H"
drought_2H_summary2_1$construct <- "Transgenic"
drought_2H_summary2_1$Day <- substr(drought_2H_summary2_1$Day,4,5)
drought_2H_summary2_1 <- as.data.frame(as.matrix(drought_2H_summary2_1[,])) 

#2
drought_2H_summary2_2 <- pivot_longer(drought_2H_summary2, cols = c(3,9,15,21), names_to = "Day",values_to = "percent_2")
drought_2H_summary2_2 <- drought_2H_summary2_2[,21:22]
drought_2H_summary2_2$event <- "2H"
drought_2H_summary2_2$construct <- "Transgenic"
drought_2H_summary2_2$Day <- substr(drought_2H_summary2_2$Day,4,5)
drought_2H_summary2_2 <- as.data.frame(as.matrix(drought_2H_summary2_2[,])) 

#3
drought_2H_summary2_3 <- pivot_longer(drought_2H_summary2, cols = c(4,10,16,22), names_to = "Day",values_to = "percent_3")
drought_2H_summary2_3 <- drought_2H_summary2_3[,21:22]
drought_2H_summary2_3$event <- "2H"
drought_2H_summary2_3$construct <- "Transgenic"
drought_2H_summary2_3$Day <- substr(drought_2H_summary2_3$Day,4,5)
drought_2H_summary2_3 <- as.data.frame(as.matrix(drought_2H_summary2_3[,])) 
#4
drought_2H_summary2_4 <- pivot_longer(drought_2H_summary2, cols = c(5,11,17,23), names_to = "Day",values_to = "percent_4")
drought_2H_summary2_4 <- drought_2H_summary2_4[,21:22]
drought_2H_summary2_4$event <- "2H"
drought_2H_summary2_4$construct <- "Transgenic"
drought_2H_summary2_4$Day <- substr(drought_2H_summary2_4$Day,4,5)
drought_2H_summary2_4 <- as.data.frame(as.matrix(drought_2H_summary2_4[,])) 

#5
drought_2H_summary2_5 <- pivot_longer(drought_2H_summary2, cols = c(6,11,18,24), names_to = "Day",values_to = "percent_5")
drought_2H_summary2_5 <- drought_2H_summary2_5[,21:22]
drought_2H_summary2_5$event <- "2H"
drought_2H_summary2_5$construct <- "Control"
drought_2H_summary2_5$Day <- substr(drought_2H_summary2_5$Day,4,5)
drought_2H_summary2_5 <- as.data.frame(as.matrix(drought_2H_summary2_5[,])) 

drought_2H_summary_3 <- inner_join(drought_2H_summary2_0, drought_2H_summary2_1, by = c("Day","event","construct"))
drought_2H_summary_3 <- inner_join(drought_2H_summary_3, drought_2H_summary2_2, by = c("Day","event","construct"))
drought_2H_summary_3 <- inner_join(drought_2H_summary_3, drought_2H_summary2_3, by = c("Day","event","construct"))
drought_2H_summary_3 <- inner_join(drought_2H_summary_3, drought_2H_summary2_4, by = c("Day","event","construct"))
#drought_2H_summary_3 <- inner_join(drought_2H_summary_3, drought_2H_summary2_5, by = c("Day","event","construct"))
drought_2H_summary_3$percent_5 <- drought_2H_summary2_5$percent_5

#####################################################################################
#### Combine event summaries together ###
Drought_summary <- rbind(drought_CT3_summary_3,drought_16_20_summary_3, drought_8_9D_summary_3, drought_5A_summary_3, drought_4A_summary_3, drought_13_15E_summary_3, drought_2H_summary_3)

write.csv(file = "LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/drought_summary_long.csv", Drought_summary)
