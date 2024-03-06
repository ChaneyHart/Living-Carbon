#script for getting daily weather t data from LC Atmos station, 2023
#install.packages("janitor")
library(janitor)
library(dplyr)
library(xts)
library(tsbox)
library(forecast)
library(lubridate)

#read in data
weather_2023 <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_weather.csv", header = TRUE, stringsAsFactors = FALSE)
weather_2022 <- read.csv("LC_2023/2023_weather/weather_processed/LC_2022_weather.csv", header = TRUE, stringsAsFactors = FALSE)
#fix timestamp format
weather_2023$Datetime <- ifelse(nchar(weather_2023$Datetime) == 10, paste0(weather_2023$Datetime," 00:00:00"),weather_2023$Datetime)
weather_2023$Datetime <- ymd_hms(weather_2023$Datetime)
weather_2022$Datetime <- ifelse(nchar(weather_2022$Datetime) == 10, paste0(weather_2022$Datetime," 00:00:00"),weather_2022$Datetime)
weather_2022$Datetime <- ymd_hms(weather_2022$Datetime)

#2022
temp_set_2022 <- subset(weather_2022, select = c("Datetime","air_temp"))
max_temp_2022 <- apply.daily(temp_set_2022, max)
min_temp_2022 <- apply.daily(temp_set_2022, min)
mean_temp_2022 <- apply.daily(temp_set_2022, mean)

VPD_set_2022 <- subset(weather_2022, select = c("Datetime","VPD"))
max_VPD_2022 <- apply.daily(VPD_set_2022, max)
min_VPD_2022 <- apply.daily(VPD_set_2022, min)


precip_set_2022 <- subset(weather_2022, select = c("Datetime","precip"))
daily_precip_2022 <- apply.daily(precip_set_2022, sum)

weather_daily_2022 <- cbind(max_temp_2022,min_temp_2022,mean_temp_2022,max_VPD_2022,min_VPD_2022,daily_precip_2022)
colnames(weather_daily_2022) <- c("max_temp","min_temp","mean_temp","max_VPD","min_VPD","daily_precip")
d2022 <- weather_daily_2022
Datetime2022 <- rownames(d2022)
rownames(d2022) <- NULL
weather_daily_2022 <- cbind(Datetime2022,d2022)
weather_daily_2022 <- weather_daily_2022 %>% mutate(Datetime = ymd_hms(Datetime2022))
weather_daily_2022$Date <- date(weather_daily_2022$Datetime)
str(weather_daily_2022)
weather_daily_2022 <- weather_daily_2022[,c(9,2,3,4,5,6,7)]
write.csv(weather_daily_2022, "LC_2023/2023_weather/weather_processed/weather_daily_2022.csv")

##read in prism data for 2022
library(stringi)
prism_2022 <- read.csv(file = "LC_2023/2023_weather/weather_raw/PRISM_dat_2022.csv",skip = 10)
prism_2022$Date <- mdy(prism_2022$Date)
str(prism_2022)
colnames(prism_2022) <- c("Date","daily_precip","min_temp","mean_temp","max_temp","min_VPD","max_VPD")
prism_2022$min_VPD <- prism_2022$min_VPD/10
prism_2022$max_VPD <- prism_2022$max_VPD/10
prism_2022 <- prism_2022[,c(1,5,3,4,7,6,2)]
#replace last day (duplicate)
prism_2022 <- prism_2022[-70,]
full_2022 <- rbind(prism_2022,weather_daily_2022)
write.csv(full_2022, file = "LC_2023/2023_weather/weather_processed/daily_weather_2022_plus_prism.csv")

#read in PRISM data to imput 2023 precip

PRISM_2023 <- read.csv(file = "LC_2023/2023_weather/PRISM_LC_early_2023.csv",skip=10)
PRISM_2023$Date <- ymd(PRISM_2023$Date)
colnames(PRISM_2023) <- c("Date","daily_precip","min_temp","mean_temp","max_temp","min_VPD","max_VPD")
PRISM_2023_precip <- subset(PRISM_2023, select = c("Date","daily_precip"))

#2023
temp_set_2023 <- subset(weather_2023, select = c("Datetime","air_temp"))
max_temp_2023 <- apply.daily(temp_set_2023, max)
min_temp_2023 <- apply.daily(temp_set_2023, min)
mean_temp_2023 <- apply.daily(temp_set_2023, mean)

VPD_set_2023 <- subset(weather_2023, select = c("Datetime","VPD"))
max_VPD_2023 <- apply.daily(VPD_set_2023, max)
min_VPD_2023 <- apply.daily(VPD_set_2023, min)


precip_set_2023 <- subset(weather_2023, select = c("Datetime","precip"))
daily_precip_2023 <- apply.daily(precip_set_2023, sum)


weather_daily_2023 <- cbind(max_temp_2023,min_temp_2023,mean_temp_2023,max_VPD_2023,min_VPD_2023,daily_precip_2023)
colnames(weather_daily_2023) <- c("max_temp","min_temp","mean_temp_2023","max_VPD","min_VPD","daily_precip")
d2023 <- weather_daily_2023
Datetime2023 <- rownames(d2023)
rownames(d2023) <- NULL
weather_daily_2023 <- cbind(Datetime2023,d2023)
weather_daily_2023 <- weather_daily_2023 %>% mutate(Datetime = ymd_hms(Datetime2023))
weather_daily_2023$Date <- as.Date(weather_daily_2023$Datetime)

#imput PRISM precip data
weather_daily_2023[1:10,7] <- PRISM_2023_precip[11:20,2]
weather_daily_2023[11:44,7] <- PRISM_2023_precip[28:61,2]


write.csv(weather_daily_2023, "LC_2023/2023_weather/weather_processed/weather_daily_2023.csv")



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
daily_predawn$Date <- as.Date(daily_predawn$Datetime)
write.csv(daily_predawn, file = "LC_2023/2023_weather/weather_processed/daily_predawn_WP.csv")


#Daily soil moisture
SM_dat <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_SM.csv")
SM_dat <- SM_dat[,c(2,3)]
str(SM_dat)
SM_dat$Datetime <- ymd_hms(SM_dat$Datetime)
SM_daily <- apply.daily(SM_dat, mean)


SM_2023 <- SM_daily
SMd2023 <- rownames(SM_daily)
rownames(SMd2023) <- NULL
SM_daily_2023 <- cbind(SMd2023,SM_2023)
SM_daily_2023$Date <- as.Date(SM_daily_2023$SMd2023)

write.csv(SM_daily_2023, file = "LC_2023/2023_weather/weather_processed/daily_SM.csv")
