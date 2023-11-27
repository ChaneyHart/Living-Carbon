#2023 weather graphing

library(ggplot2)
library(lubridate)

#read in data
weather_tot <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_weather.csv")
weather_tot$Datetime <- ymd_hms(weather_tot$Datetime)
weather_hourly <- read.csv("LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
weather_hourly$Datetime <- ymd_hms(weather_hourly$Datetime)
weather_daily <- read.csv("LC_2023/2023_weather/weather_processed/weather_daily_2023.csv")
weather_daily$Datetime <- ymd_hms(weather_daily$Datetime)
daily_predawn <- read.csv("LC_2023/2023_weather/weather_processed/daily_predawn_WP.csv")
daily_predawn$Datetime <- ymd_hms(daily_predawn$Datetime)
daily_predawn$psi_pd <- (-1*daily_predawn$psi_pd)
SM_dat <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_SM.csv")
SM_dat$Datetime <- ymd_hms(SM_dat$Datetime)

#graph temp
str(weather_tot)
temp_plot <- ggplot(weather_daily, aes(x=Datetime))+
  geom_line(aes(y=max_temp, color = "red4"))+
  geom_line(aes(y=min_temp),color="pink")+
  ylab("Daily max and min air temperature ˚C")+
  theme_bw()
temp_plot
ggsave(filename = "LC_2023/2023_weather/weather_graphs/temp_plot.png", plot = temp_plot, width = 8,height = 2,units = "in",dpi = 300)

#graph VPD
VPD_plot <- ggplot(weather_daily, aes(x=Datetime))+
  geom_line(aes(y=max_VPD),color="black")+
  geom_line(aes(y=min_VPD),color="gray")+
  ylab("Daily max and min VPD (KPa)")+
  theme_bw()

VPD_plot
ggsave(filename = "LC_2023/2023_weather/weather_graphs/VPD_plot.png", plot = VPD_plot, width = 8,height = 2,units = "in",dpi = 300)
#graph precip

precip_plot<- ggplot(weather_daily, aes(x=Datetime))+
  geom_line(aes(y=daily_precip),color="blue4")+
  ylab("Daily Precipiation (mm)")+
  theme_bw()
precip_plot

ggsave(filename = "LC_2023/2023_weather/weather_graphs/precip_plot.png", plot = precip_plot, width = 8,height = 2,units = "in",dpi = 300)
#graph SM and WP


Soil_plot <- ggplot()+
  geom_line(data = subset(SM_dat, Datetime > as.POSIXct("2023-07-20 00:00") & Datetime < as.POSIXct("2023-10-11 00:00")),aes(x= Datetime, y=soil_water_content),color="cyan3")+
  geom_point(data = daily_predawn, aes(x=Datetime, y=psi_pd/10),color="black")+
  geom_line(data = daily_predawn, aes(x=Datetime, y=psi_pd/10),color="black")+
  scale_y_continuous(
    # Features of the first axis
    name = "30 cm Soil water content (m3/m3)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*10, name="Est. soil water potential (KPa)")
  ) +
  theme_bw()

Soil_plot

ggsave(filename = "LC_2023/2023_weather/weather_graphs/soil_plot.png", plot = Soil_plot, width = 5,height = 2,units = "in",dpi = 300)

