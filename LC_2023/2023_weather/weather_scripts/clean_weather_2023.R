#script to match weather station data with A-Ci measurements

#read in csv with weather readings for growing season of 2023
weather_dat_2023 <- read.csv("")

#find averages for unique date time combos
air_temp <- as.data.frame(aggregate(X.C.Air.Temperature ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
VPD <- as.data.frame(aggregate(kPa.Vapor.Pressure ~ mdh, data = weather_dat, FUN=mean, na.rm=TRUE))
