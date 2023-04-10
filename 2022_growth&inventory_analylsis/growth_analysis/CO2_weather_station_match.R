#script to match weather station data with A-Ci measurements

#read in csv with weather readings for the month of July 2022
weather_dat <- read.csv("2022_physiology_analysis/Weather_station/Aci_Atmos41_data.csv")

#find averages for unique date time combos
air_temp <- as.data.frame(aggregate(X.C.Air.Temperature ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
VPD <- as.data.frame(aggregate(kPa.Vapor.Pressure ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
