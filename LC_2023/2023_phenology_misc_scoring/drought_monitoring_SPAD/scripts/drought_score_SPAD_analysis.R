#LC growth analysis update - Sep 21 2023


library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)

#read in data

drought_scores_SPAD <- read.csv("LC_2023/2023_phenology_misc_scoring/Drought_Scores_SPAD_compiled.csv")
growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
drought_scores_SPAD <- inner_join(drought_scores_SPAD,growth)

#seperate so NAs in one do not affect the other
drought_scores_2023 <- subset(drought_scores_SPAD, select = c("row","column","event_short","construct","construct2","ID","DS_1","DS_18","DS_32","DS_55"))
SPAD_2023 <- subset(drought_scores_SPAD, select = c("row","column","event_short","construct","construct2","ID","SPAD_18","SPAD_32","SPAD_55"))

###########################################
#### cleaning drought scores
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


###########################
### summarize proproprtions of scores for each event

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


drought_13_15E <- subset(drought_scores_2023, event_short == "13-15E")

drought_13_15E_summary <- drought_13_15E %>% summarize(
  day1_0 = count(subset(drought_13_15E, DS_1 == "0")),
  day1_1 = count(subset(drought_13_15E, DS_1 == "1")),
  day1_2 = count(subset(drought_13_15E, DS_1 == "2")),
  day1_3 = count(subset(drought_13_15E, DS_1 == "3")),
  day1_4 = count(subset(drought_13_15E, DS_1 == "4")),
  day1_5 = count(subset(drought_13_15E, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_13_15E, DS_18 == "0")),
  day18_1 = count(subset(drought_13_15E, DS_18 == "1")),
  day18_2 = count(subset(drought_13_15E, DS_18 == "2")),
  day18_3 = count(subset(drought_13_15E, DS_18 == "3")),
  day18_4 = count(subset(drought_13_15E, DS_18 == "4")),
  day18_5 = count(subset(drought_13_15E, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_13_15E, DS_32 == "0")),
  day32_1 = count(subset(drought_13_15E, DS_32 == "1")),
  day32_2 = count(subset(drought_13_15E, DS_32 == "2")),
  day32_3 = count(subset(drought_13_15E, DS_32 == "3")),
  day32_4 = count(subset(drought_13_15E, DS_32 == "4")),
  day32_5 = count(subset(drought_13_15E, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_13_15E, DS_55 == "0")),
  day55_1 = count(subset(drought_13_15E, DS_55 == "1")),
  day55_2 = count(subset(drought_13_15E, DS_55 == "2")),
  day55_3 = count(subset(drought_13_15E, DS_55 == "3")),
  day55_4 = count(subset(drought_13_15E, DS_55 == "4")),
  day55_5 = count(subset(drought_13_15E, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_13_15E_summary2 <- drought_13_15E_summary[1,c(7,14,21,28)]
drought_13_15E_summary2$event <- "13-15E"
drought_13_15E_summary2$construct <- "transgenic"
as.data.frame(as.matrix(drought_13_15E_summary2[c(1,2,3,4,5,6)])) 
drought_13_15E_summary2 <- as.data.frame(as.matrix(drought_13_15E_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_13_15E_summary2) <- c("1","18","32","55","event","construct")


drought_2H <- subset(drought_scores_2023, event_short == "2H")

drought_2H_summary <- drought_2H %>% summarize(
  day1_0 = count(subset(drought_2H, DS_1 == "0")),
  day1_1 = count(subset(drought_2H, DS_1 == "1")),
  day1_2 = count(subset(drought_2H, DS_1 == "2")),
  day1_3 = count(subset(drought_2H, DS_1 == "3")),
  day1_4 = count(subset(drought_2H, DS_1 == "4")),
  day1_5 = count(subset(drought_2H, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_2H, DS_18 == "0")),
  day18_1 = count(subset(drought_2H, DS_18 == "1")),
  day18_2 = count(subset(drought_2H, DS_18 == "2")),
  day18_3 = count(subset(drought_2H, DS_18 == "3")),
  day18_4 = count(subset(drought_2H, DS_18 == "4")),
  day18_5 = count(subset(drought_2H, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_2H, DS_32 == "0")),
  day32_1 = count(subset(drought_2H, DS_32 == "1")),
  day32_2 = count(subset(drought_2H, DS_32 == "2")),
  day32_3 = count(subset(drought_2H, DS_32 == "3")),
  day32_4 = count(subset(drought_2H, DS_32 == "4")),
  day32_5 = count(subset(drought_2H, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_2H, DS_55 == "0")),
  day55_1 = count(subset(drought_2H, DS_55 == "1")),
  day55_2 = count(subset(drought_2H, DS_55 == "2")),
  day55_3 = count(subset(drought_2H, DS_55 == "3")),
  day55_4 = count(subset(drought_2H, DS_55 == "4")),
  day55_5 = count(subset(drought_2H, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_2H_summary2 <- drought_2H_summary[1,c(7,14,21,28)]
drought_2H_summary2$event <- "2H"
drought_2H_summary2$construct <- "transgenic"
as.data.frame(as.matrix(drought_2H_summary2[c(1,2,3,4,5,6)])) 
drought_2H_summary2 <- as.data.frame(as.matrix(drought_2H_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_2H_summary2) <- c("1","18","32","55","event","construct")


drought_5A <- subset(drought_scores_2023, event_short == "5A")

drought_5A_summary <- drought_5A %>% summarize(
  day1_0 = count(subset(drought_5A, DS_1 == "0")),
  day1_1 = count(subset(drought_5A, DS_1 == "1")),
  day1_2 = count(subset(drought_5A, DS_1 == "2")),
  day1_3 = count(subset(drought_5A, DS_1 == "3")),
  day1_4 = count(subset(drought_5A, DS_1 == "4")),
  day1_5 = count(subset(drought_5A, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_5A, DS_18 == "0")),
  day18_1 = count(subset(drought_5A, DS_18 == "1")),
  day18_2 = count(subset(drought_5A, DS_18 == "2")),
  day18_3 = count(subset(drought_5A, DS_18 == "3")),
  day18_4 = count(subset(drought_5A, DS_18 == "4")),
  day18_5 = count(subset(drought_5A, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_5A, DS_32 == "0")),
  day32_1 = count(subset(drought_5A, DS_32 == "1")),
  day32_2 = count(subset(drought_5A, DS_32 == "2")),
  day32_3 = count(subset(drought_5A, DS_32 == "3")),
  day32_4 = count(subset(drought_5A, DS_32 == "4")),
  day32_5 = count(subset(drought_5A, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_5A, DS_55 == "0")),
  day55_1 = count(subset(drought_5A, DS_55 == "1")),
  day55_2 = count(subset(drought_5A, DS_55 == "2")),
  day55_3 = count(subset(drought_5A, DS_55 == "3")),
  day55_4 = count(subset(drought_5A, DS_55 == "4")),
  day55_5 = count(subset(drought_5A, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_5A_summary2 <- drought_5A_summary[1,c(7,14,21,28)]
drought_5A_summary2$event <- "5A"
drought_5A_summary2$construct <- "transgenic"
drought_5A_summary2 <- as.data.frame(as.matrix(drought_5A_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_5A_summary2) <- c("1","18","32","55","event","construct")

drought_4A <- subset(drought_scores_2023, event_short == "4A")

drought_4A_summary <- drought_4A %>% summarize(
  day1_0 = count(subset(drought_4A, DS_1 == "0")),
  day1_1 = count(subset(drought_4A, DS_1 == "1")),
  day1_2 = count(subset(drought_4A, DS_1 == "2")),
  day1_3 = count(subset(drought_4A, DS_1 == "3")),
  day1_4 = count(subset(drought_4A, DS_1 == "4")),
  day1_5 = count(subset(drought_4A, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_4A, DS_18 == "0")),
  day18_1 = count(subset(drought_4A, DS_18 == "1")),
  day18_2 = count(subset(drought_4A, DS_18 == "2")),
  day18_3 = count(subset(drought_4A, DS_18 == "3")),
  day18_4 = count(subset(drought_4A, DS_18 == "4")),
  day18_5 = count(subset(drought_4A, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_4A, DS_32 == "0")),
  day32_1 = count(subset(drought_4A, DS_32 == "1")),
  day32_2 = count(subset(drought_4A, DS_32 == "2")),
  day32_3 = count(subset(drought_4A, DS_32 == "3")),
  day32_4 = count(subset(drought_4A, DS_32 == "4")),
  day32_5 = count(subset(drought_4A, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_4A, DS_55 == "0")),
  day55_1 = count(subset(drought_4A, DS_55 == "1")),
  day55_2 = count(subset(drought_4A, DS_55 == "2")),
  day55_3 = count(subset(drought_4A, DS_55 == "3")),
  day55_4 = count(subset(drought_4A, DS_55 == "4")),
  day55_5 = count(subset(drought_4A, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_4A_summary2 <- drought_4A_summary[1,c(7,14,21,28)]
drought_4A_summary2$event <- "4A"
drought_4A_summary2$construct <- "transgenic"
drought_4A_summary2 <- as.data.frame(as.matrix(drought_4A_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_4A_summary2) <- c("1","18","32","55","event","construct")


drought_16_20 <- subset(drought_scores_2023, event_short == "16-20")

drought_16_20_summary <- drought_16_20 %>% summarize(
  day1_0 = count(subset(drought_16_20, DS_1 == "0")),
  day1_1 = count(subset(drought_16_20, DS_1 == "1")),
  day1_2 = count(subset(drought_16_20, DS_1 == "2")),
  day1_3 = count(subset(drought_16_20, DS_1 == "3")),
  day1_4 = count(subset(drought_16_20, DS_1 == "4")),
  day1_5 = count(subset(drought_16_20, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_16_20, DS_18 == "0")),
  day18_1 = count(subset(drought_16_20, DS_18 == "1")),
  day18_2 = count(subset(drought_16_20, DS_18 == "2")),
  day18_3 = count(subset(drought_16_20, DS_18 == "3")),
  day18_4 = count(subset(drought_16_20, DS_18 == "4")),
  day18_5 = count(subset(drought_16_20, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_16_20, DS_32 == "0")),
  day32_1 = count(subset(drought_16_20, DS_32 == "1")),
  day32_2 = count(subset(drought_16_20, DS_32 == "2")),
  day32_3 = count(subset(drought_16_20, DS_32 == "3")),
  day32_4 = count(subset(drought_16_20, DS_32 == "4")),
  day32_5 = count(subset(drought_16_20, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_16_20, DS_55 == "0")),
  day55_1 = count(subset(drought_16_20, DS_55 == "1")),
  day55_2 = count(subset(drought_16_20, DS_55 == "2")),
  day55_3 = count(subset(drought_16_20, DS_55 == "3")),
  day55_4 = count(subset(drought_16_20, DS_55 == "4")),
  day55_5 = count(subset(drought_16_20, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_16_20_summary2 <- drought_16_20_summary[1,c(7,14,21,28)]
drought_16_20_summary2$event <- "16-20"
drought_16_20_summary2$construct <- "escape"
drought_16_20_summary2 <- as.data.frame(as.matrix(drought_16_20_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_16_20_summary2) <- c("1","18","32","55","event","construct")


drought_8_9D <- subset(drought_scores_2023, event_short == "8-9D")

drought_8_9D_summary <- drought_8_9D %>% summarize(
  day1_0 = count(subset(drought_8_9D, DS_1 == "0")),
  day1_1 = count(subset(drought_8_9D, DS_1 == "1")),
  day1_2 = count(subset(drought_8_9D, DS_1 == "2")),
  day1_3 = count(subset(drought_8_9D, DS_1 == "3")),
  day1_4 = count(subset(drought_8_9D, DS_1 == "4")),
  day1_5 = count(subset(drought_8_9D, DS_1 == "5")),
  day1_prop_0 = day1_0/(sum(day1_0,day1_1,day1_2,day1_3,day1_4,day1_5)),
  day18_0 = count(subset(drought_8_9D, DS_18 == "0")),
  day18_1 = count(subset(drought_8_9D, DS_18 == "1")),
  day18_2 = count(subset(drought_8_9D, DS_18 == "2")),
  day18_3 = count(subset(drought_8_9D, DS_18 == "3")),
  day18_4 = count(subset(drought_8_9D, DS_18 == "4")),
  day18_5 = count(subset(drought_8_9D, DS_18 == "5")),
  day18_prop_0 = day18_0/(sum(day18_0,day18_1,day18_2,day18_3,day18_4,day18_5)),
  day32_0 = count(subset(drought_8_9D, DS_32 == "0")),
  day32_1 = count(subset(drought_8_9D, DS_32 == "1")),
  day32_2 = count(subset(drought_8_9D, DS_32 == "2")),
  day32_3 = count(subset(drought_8_9D, DS_32 == "3")),
  day32_4 = count(subset(drought_8_9D, DS_32 == "4")),
  day32_5 = count(subset(drought_8_9D, DS_32 == "5")),
  day32_prop_0 = day32_0/(sum(day32_0,day32_1,day32_2,day32_3,day32_4,day32_5)),
  day55_0 = count(subset(drought_8_9D, DS_55 == "0")),
  day55_1 = count(subset(drought_8_9D, DS_55 == "1")),
  day55_2 = count(subset(drought_8_9D, DS_55 == "2")),
  day55_3 = count(subset(drought_8_9D, DS_55 == "3")),
  day55_4 = count(subset(drought_8_9D, DS_55 == "4")),
  day55_5 = count(subset(drought_8_9D, DS_55 == "5")),
  day55_prop_0 = day55_0/(sum(day55_0,day55_1,day55_2,day55_3,day55_4,day55_5))
  
)
drought_8_9D_summary2 <- drought_8_9D_summary[1,c(7,14,21,28)]
drought_8_9D_summary2$event <- "8-9D"
drought_8_9D_summary2$transgenic <- "escape"
drought_8_9D_summary2 <- as.data.frame(as.matrix(drought_8_9D_summary2[c(1,2,3,4,5,6)])) 
colnames(drought_8_9D_summary2) <- c("1","18","32","55","event","construct")


drought_summary <- bind_rows(drought_CT3_summary2,drought_13_15E_summary2,drought_2H_summary2,drought_5A_summary2,drought_4A_summary2,drought_16_20_summary2,drought_8_9D_summary2)
#phen.top_long <- pivot_longer(phen.top, cols = c(Top_0,Top_2,Top_6), names_to = "Day", values_to = "Score")
drought_summary_long <-pivot_longer(drought_summary, cols = c(1,2,3,4), names_to = "Day", values_to = "proportion")

str(drought_summary_long)
drought_summary_long$Day <- as.numeric(drought_summary_long$Day)
drought_summary_long$proportion <- as.numeric(drought_summary_long$proportion)


library(ggrepel)

drought_timeline <- ggplot(drought_summary_long, aes(x=Day,y=proportion,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of non-drought trees (recieved score 0)")+
  theme_bw()

drought_timeline2 <- drought_timeline + colorScale
drought_timeline2
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_timeline.png",plot= drought_timeline2, width = 6, height = 5, units = "in", dpi=300)

##################################################
#SPAD
SPAD_2023$SPAD_18 <- as.numeric(SPAD_2023$SPAD_18)
SPAD_2023$SPAD_32 <- as.numeric(SPAD_2023$SPAD_32)
SPAD_2023$SPAD_55 <- as.numeric(SPAD_2023$SPAD_55)
SPAD_2023 <- na.omit(SPAD_2023)

#graphs and data tables


SPAD_subset <- subset(SPAD_2023, event_short == "16-20"|event_short == "8-9D"|event_short == "5A"|event_short == "4A"|event_short == "2H"|event_short == "13-15E"| event_short == "CT3")

SPAD_long <- pivot_longer(SPAD_subset, cols = (c(7,8,9)), names_to = "Day", values_to = "SPAD")
SPAD_long <- SPAD_long %>% mutate(Day = as.numeric(substr(SPAD_long$Day,6,7)))


SPAD_summary_full <- SPAD_long %>% group_by(event_short,Day) %>% summarize(
  n = n(),
  sd_SPAD = sd(SPAD),
  mean_SPAD = mean(SPAD),
  se_SPAD = sd_SPAD/sqrt(n),
  tst_SPAD = qt(0.975,(n-1)),
  upper_SPAD = mean_SPAD + tst_SPAD*se_SPAD,
  lower_SPAD = mean_SPAD - tst_SPAD*se_SPAD
)

SPAD_summary_full <- SPAD_summary_full %>% mutate(construct = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A"|event_short == "13-15E" | event_short == "2H" ~"transgenic",
  event_short == "16-20" | event_short == "8-9D"~"escape",
  event_short == "CT3" ~ "control"))

SPAD_plot <- ggplot(SPAD_summary_full, aes(x=Day,y= mean_SPAD, color = construct, shape = event_short))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin=lower_SPAD,ymax=upper_SPAD),width =0.7,alpha =0.7,linetype="dashed")+
  ylab("chlorophyll content (SPAD)")+
  xlab("Day (since July 20th")+
  theme_bw()

SPAD_plot_2 <- SPAD_plot + colorScale

ggsave(filename = "LC_2023/2023_phenology_misc_scoring/SPAD_plot_full.png",plot=SPAD_plot_2, dpi = 300)


SPAD_summary <- data_frame()
event_SPAD <- as.data.frame(aggregate(SPAD_18 ~ event_short, data = SPAD_2023, FUN = mean))
event_var <- aggregate(SPAD_18 ~ event_short, data = SPAD_2023, FUN = var)
event_n <- aggregate(SPAD_18 ~ event_short, data = SPAD_2023, FUN = length)

event_SPAD_summary <- inner_join(event_SPAD, event_var, by = "event_short")
event_SPAD_summary <- inner_join(event_SPAD_summary, event_n, by = "event_short")

colnames(event_SPAD_summary) <- c("Event", "mean_SPAD","var","n")

event_SPAD_summary$se <- sqrt(event_SPAD_summary$var/event_SPAD_summary$n)
event_SPAD_summary$tst <- qt(0.975, event_SPAD_summary$n)
event_SPAD_summary$sd <- sqrt(event_SPAD_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_SPAD_summary$upper_SPAD <- event_SPAD_summary$mean_SPAD + (event_SPAD_summary$tst*event_SPAD_summary$se)
event_SPAD_summary$lower_SPAD <- event_SPAD_summary$mean_SPAD - (event_SPAD_summary$tst*event_SPAD_summary$se)

write.csv(event_SPAD_summary, file = "LC_2023/2023_growth_inventory_analysis/event_SPAD_table.csv")

event_SPAD_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

myColors <- c('gray',
              'cyan4',
              'red4')

names(myColors) <- levels(event_SPAD_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)
SPAD_plot<- ggplot(subset(event_SPAD_summary, Event == "16-20"|Event == "8-9D"|Event == "5A"|Event == "4A"|Event == "2H"|Event == "13-15E"| Event == "CT3"), aes(x=reorder(Event,mean_SPAD),y=mean_SPAD,color=construct))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_SPAD,ymax=upper_SPAD))+
  xlab("Event")+
  ylab("Chlorophyll content (SPAD)")

SPAD_plot + colorScale


SPAD_summary2 <- data_frame()
event_SPAD2 <- as.data.frame(aggregate(SPAD_32 ~ event_short, data = SPAD_2023, FUN = mean))
event_var2 <- aggregate(SPAD_32 ~ event_short, data = SPAD_2023, FUN = var)
event_n2 <- aggregate(SPAD_32 ~ event_short, data = SPAD_2023, FUN = length)

event_SPAD_summary2 <- inner_join(event_SPAD2, event_var, by = "event_short")
event_SPAD_summary2 <- inner_join(event_SPAD_summary2, event_n, by = "event_short")

event_SPAD_summary2$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

colnames(event_SPAD_summary2) <- c("Event", "mean_SPAD","var","n")

event_SPAD_summary2$se <- sqrt(event_SPAD_summary2$var/event_SPAD_summary2$n)
event_SPAD_summary2$tst <- qt(0.975, event_SPAD_summary$n)
event_SPAD_summary2$sd <- sqrt(event_SPAD_summary2$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_SPAD_summary2$upper_SPAD <- event_SPAD_summary2$mean_SPAD + (event_SPAD_summary2$tst*event_SPAD_summary2$se)
event_SPAD_summary2$lower_SPAD <- event_SPAD_summary2$mean_SPAD - (event_SPAD_summary2$tst*event_SPAD_summary2$se)

write.csv(event_SPAD_summary2, file = "LC_2023/2023_growth_inventory_analysis/event_SPAD_table.csv")

event_SPAD_summary2$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

myColors <- c('gray',
              'cyan4',
              'red4')


names(myColors) <- levels(event_SPAD_summary2$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)

SPAD_plot2<- ggplot(subset(event_SPAD_summary2,Event == "16-20"|Event == "8-9D"|Event == "5A"|Event == "4A"|Event == "2H"|Event == "13-15E"| Event == "CT3"), aes(x=reorder(Event,mean_SPAD),y=mean_SPAD,color=construct))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_SPAD,ymax=upper_SPAD))+
  xlab("Event")+
  ylab("SPAD")

SPAD_plot2 + colorScale


