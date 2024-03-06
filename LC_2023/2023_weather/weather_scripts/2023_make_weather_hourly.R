#script for getting average hourly temperature and VPD data from LC Atmos station, 2023
#install.packages("janitor")
library(janitor)
library(dplyr)

#read in data
weather <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_weather.csv", header = TRUE, stringsAsFactors = FALSE)

weather_subset <- subset(weather, select =c(Datetime,air_temp, VPD,solar_radiation))

weather_subset$Datetime <- ymd_hms(weather_subset$Datetime)
str(weather_subset)

hourseq = seq.POSIXt(min(weather_subset$Datetime), max(weather_subset$Datetime), by='1 hour')
weather_hourly = weather_subset %>% group_by(Hourly = cut(Datetime, breaks=hourseq)) %>%
  summarise_each(list(~mean(., na.rm=T), ~max(., na.rm=T)))
weather_hourly <- subset(weather_hourly, select = c("Hourly","air_temp_mean","VPD_mean","air_temp_max","VPD_max","solar_radiation_mean","solar_radiation_max"))

colnames(weather_hourly) <- c("Datetime","air_temp_mean","VPD_mean","air_temp_max","VPD_max","solar_radiation_mean","solar_radiation_max")

write.csv(weather_hourly, "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")



