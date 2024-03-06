#script to match weather station data with A-Ci measurements
library(dplyr)
library(xts)
library(tsbox)
library(forecast)
library(lubridate)

#read in csv with weather readings for growing season of 2023
weather_dat_2023_1 <- read.csv("LC_2023/2023_weather/weather_raw/atmos_2023_growing_Season_4_16_7_20.csv",skip = 2, header = T,sep = ',', stringsAsFactors = F)
weather_dat_2023_2 <- read.csv("LC_2023/2023_weather/weather_raw/atmos_2023_growing_Season_7_20_10_16.csv",skip = 2, header = T,sep = ',', stringsAsFactors = F)

#read in csv with weather dat for 2022

weather_dat_2022 <- read.csv("LC_2023/2023_weather/weather_raw/atmos_2022_gro_seaspm.csv",skip = 2, header = T,sep = ',', stringsAsFactors = F)

#modify column names
colnames(weather_dat_2023_1) <- c("Datetime","solar_radiation","precip","lightning","lightning_dist","wind_dir","wind_speed","gust_speed","air_temp","RH","atmos_pressure","X_level","Y_level","max_precip_rate","RH_temp","batt_percent","batt_volt","ref_pressure","logger_temp")
colnames(weather_dat_2023_2) <- c("Datetime","solar_radiation","precip","lightning","lightning_dist","wind_dir","wind_speed","gust_speed","air_temp","RH","atmos_pressure","X_level","Y_level","max_precip_rate","RH_temp","soil_water_content","batt_percent","batt_volt","ref_pressure","logger_temp")
colnames(weather_dat_2022) <- c("Datetime","solar_radiation","precip","lightning","lightning_dist","wind_dir","wind_speed","gust_speed","air_temp","RH","atmos_pressure","X_level","Y_level","max_precip_rate","RH_temp","batt_percent","batt_volt","ref_pressure","logger_temp")

#convert date time to correct format
weather_dat_2023_1 <- weather_dat_2023_1 %>% mutate(Datetime = mdy_hm(Datetime))
weather_dat_2023_2 <- weather_dat_2023_2 %>% mutate(Datetime = mdy_hm(Datetime))
weather_dat_2022 <- weather_dat_2022 %>% mutate(Datetime = mdy_hm(Datetime))

#create seperate contiuous SM and weather datasets
weather_dat_2023_2_SM <- subset(weather_dat_2023_2, select = c("Datetime","soil_water_content"))
weather_dat_2023_2 <- subset(weather_dat_2023_2, select = c("Datetime","solar_radiation","precip","lightning","lightning_dist","wind_dir","wind_speed","gust_speed","air_temp","RH","atmos_pressure","X_level","Y_level","max_precip_rate","RH_temp","batt_percent","batt_volt","ref_pressure","logger_temp"))

#join discontinuous weather meas

LC_2023_weather <- rbind(weather_dat_2023_1, weather_dat_2023_2)

str(LC_2023_weather)


#calculate VPDleaf from dataset
LC_2023_weather$es <- 0.6108*(exp((17.2694*LC_2023_weather$air_temp)/(LC_2023_weather$air_temp+237.3)))
#Tetens equation from Oregon State Dewpoint and VPD calculations pdf
#SVP is in millibars
LC_2023_weather$ea <- (LC_2023_weather$RH*LC_2023_weather$es)
LC_2023_weather$VPD <- LC_2023_weather$es-LC_2023_weather$ea

weather_dat_2022$es <- 0.6108*(exp((17.2694*weather_dat_2022$air_temp)/(weather_dat_2022$air_temp+237.3)))
#Tetens equation from Oregon State Dewpoint and VPD calculations pdf
#SVP is in millibars
weather_dat_2022$ea <- (weather_dat_2022$RH*weather_dat_2022$es)
weather_dat_2022$VPD <- weather_dat_2022$es-weather_dat_2022$ea



#output_files
write.csv(file ="LC_2023/2023_weather/LC_2023_weather.csv",LC_2023_weather)
write.csv(file ="LC_2023/2023_weather/LC_2023_SM.csv",weather_dat_2023_2_SM)
write.csv(file = "LC_2023/2023_weather/weather_processed/LC_2022_weather.csv",weather_dat_2022)
#air_temp <- as.data.frame(aggregate(X.C.Air.Temperature ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
#VPD <- as.data.frame(aggregate(kPa.Vapor.Pressure ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
