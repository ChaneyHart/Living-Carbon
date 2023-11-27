#script for getting daily weather t data from LC Atmos station, 2023
#install.packages("janitor")
library(janitor)
library(dplyr)
library(xts)
library(tsbox)
library(forecast)
library(lubridate)

#read in data
weather <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_weather.csv", header = TRUE, stringsAsFactors = FALSE)
weather$Datetime <- ymd_hms(weather$Datetime)


temp_set <- subset(weather, select = c("Datetime","air_temp"))
max_temp <- apply.daily(temp_set, max)
min_temp <- apply.daily(temp_set, min)
mean_temp <- apply.daily(temp_set, mean)

VPD_set <- subset(weather, select = c("Datetime","VPD"))
max_VPD <- apply.daily(VPD_set, max)
min_VPD <- apply.daily(VPD_set, min)


precip_set <- subset(weather, select = c("Datetime","precip"))
daily_precip <- apply.daily(precip_set, sum)

weather_daily_2023 <- cbind(max_temp,min_temp,max_VPD,min_VPD,daily_precip)
colnames(weather_daily_2023) <- c("max_temp","min_temp","max_VPD","min_VPD","daily_precip")
d <- weather_daily_2023
Datetime <- rownames(d)
rownames(d) <- NULL
weather_daily_2023 <- cbind(Datetime,d)
weather_daily_2023 <- weather_daily_2023 %>% mutate(Datetime = ymd_hms(Datetime))


write.csv(weather_daily_2023, "LC_2023/2023_weather/weather_daily_2023.csv")

#predawn daily water potentials
#from scattered days

WP_8_1 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_field_conditions/8_1_23_predawn_WP_filled.csv")
WP_8_1 <- subset(WP_8_1, select = "Predawn_MPA")
?aggregate
WP_8_1 <- WP_8_1 %>% summarize(
  psi_pd = mean(Predawn_MPA)
)
WP_8_1$Datetime <- 	ymd_hms("2023-08-01 23:45:00")

WP_8_16 <- read.csv(file = "LC_2023/2023_molecular_analysis/August_molecular/Metabolite_datasheet_8_16_conditions_filled.csv")
WP_8_16 <- subset(WP_8_16, select = "Mpa_predawn")
WP_8_16$Mpa_predawn <- as.numeric(WP_8_16$Mpa_predawn)
WP_8_16 <- na.omit(WP_8_16)
WP_8_16 <- WP_8_16 %>% summarize(
  psi_pd = mean(Mpa_predawn)
)
WP_8_16$Datetime <- 	ymd_hms("2023-08-16 23:45:00")

WP_8_17 <- read.csv(file = "LC_2023/2023_molecular_analysis/August_molecular/Metabolite_datasheet_8_17_conditions_filled.csv")
WP_8_17 <- subset(WP_8_17, select = "Mpa_predawn")
WP_8_17$Mpa_predawn <- as.numeric(WP_8_17$Mpa_predawn)
WP_8_17 <- na.omit(WP_8_17)
WP_8_17 <- WP_8_17 %>% summarize(
  psi_pd = mean(Mpa_predawn)
)
WP_8_17$Datetime <- 	ymd_hms("2023-08-17 23:45:00")


WP_8_28 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_field_conditions/August_diurnal_mPA.csv")
#exclude irrigated trees
Irrigation_set <- read.csv("LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
WP_8_28 <- WP_8_28[!(WP_8_28$ID %in% Irrigation_set$ID),]

WP_8_28 <- subset(WP_8_28, select = "mPA")
WP_8_28$mPA <- as.numeric(WP_8_28$mPA)
WP_8_28 <- na.omit(WP_8_28)
WP_8_28 <- WP_8_28 %>% summarize(
  psi_pd = mean(mPA)
)
WP_8_28$Datetime <- 	ymd_hms("2023-08-28 23:45:00")

daily_predawn <- rbind(WP_8_1,WP_8_16,WP_8_17,WP_8_28)
write.csv(daily_predawn, file = "LC_2023/2023_weather/weather_processed/daily_predawn_WP.csv")
