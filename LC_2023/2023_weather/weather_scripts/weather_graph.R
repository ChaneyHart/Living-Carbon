#2023 weather graphing

library(ggplot2)
library(lubridate)

#read in data
weather_tot <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_weather.csv")
weather_tot$Datetime <- ymd_hms(weather_tot$Datetime)
weather_hourly <- read.csv("LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
weather_hourly$Datetime <- ymd_hms(weather_hourly$Datetime)
weather_daily <- read.csv("LC_2023/2023_weather/weather_processed/weather_daily_2023.csv")
weather_daily$Date <- date(weather_daily$Datetime)
weather_daily_2022 <- read.csv("LC_2023/2023_weather/weather_processed/daily_weather_2022_plus_prism.csv")
weather_daily_2022$Date <- date(weather_daily_2022$Date)
daily_predawn <- read.csv("LC_2023/2023_weather/weather_processed/daily_predawn_WP.csv")
daily_predawn$Datetime <- ymd_hms(daily_predawn$Datetime)
daily_predawn$psi_pd <- (-1*daily_predawn$psi_pd)
SM_dat <- read.csv("LC_2023/2023_weather/weather_processed/LC_2023_SM.csv")
SM_dat$Datetime <- ymd_hms(SM_dat$Datetime)

#graph temp
str(weather_tot)
temp_plot_23 <- ggplot(weather_daily, aes(x=Date))+
  geom_line(aes(y=max_temp, color = "red4"))+
  geom_line(aes(y=min_temp),color="pink")+
  ylab("Daily max and min air temperature ˚C")+
  ggtitle(label = "2023 growing season")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

temp_plot_23
ggsave(filename = "LC_2023/2023_weather/weather_graphs/temp_plot.png", plot = temp_plot_23, width = 8,height = 2,units = "in",dpi = 300)

temp_plot_22 <- ggplot(subset(weather_daily_2022, Date < "2022-10-15"), aes(x=Date))+
  geom_line(aes(y=max_temp, color = "red4"))+
  geom_line(aes(y=min_temp),color="pink")+
  ylab("Daily max and min air temperature ˚C")+
  ggtitle(label = "2022 growing season")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

temp_plot_22
ggsave(filename = "LC_2023/2023_weather/weather_graphs/temp_plot_2022.png", plot = temp_plot_22, width = 8,height = 2,units = "in",dpi = 300)

#graph VPD
str(weather_daily)
VPD_plot <- ggplot(weather_daily, aes(x=Date))+
  geom_line(aes(y=max_VPD),color="black")+
  geom_line(aes(y=min_VPD),color="gray")+
  ylab("Daily max and min VPD (KPa)")+
  ggtitle(label = "2023 growing season")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

VPD_plot
ggsave(filename = "LC_2023/2023_weather/weather_graphs/VPD_plot.png", plot = VPD_plot, width = 8,height = 2,units = "in",dpi = 300)

VPD_plot_22 <- ggplot(subset(weather_daily_2022, Date < "2022-10-15"), aes(x=Date))+
  geom_line(aes(y=max_VPD),color="black")+
  geom_line(aes(y=min_VPD),color="gray")+
  ylab("Daily max and min VPD (KPa)")+
  theme_bw()+
  ggtitle(label = "2022 growing season")+
  theme(plot.title = element_text(hjust = 0.5))

VPD_plot_22
#graph precip

precip_plot <- ggplot(weather_daily, aes(x=Date))+
  geom_line(aes(y=daily_precip),color="blue4")+
  ylab("Daily Precipiation (mm)")+
  theme_bw()+
  ggtitle(label = "2023 growing season")+
  theme(plot.title = element_text(hjust = 0.5))
precip_plot

#need to imput weather from PRISM here
ggsave(filename = "LC_2023/2023_weather/weather_graphs/precip_plot.png", plot = precip_plot, width = 8,height = 2,units = "in",dpi = 300)
#graph SM and WP

precip_plot_2022 <- ggplot(subset(weather_daily_2022, Date < "2022-10-15"), aes(x=Date))+
  geom_line(aes(y=daily_precip),color="blue4")+
  ylab("Daily Precipiation (mm)")+
  theme_bw()+
  ggtitle(label = "2022 growing season")+
  theme(plot.title = element_text(hjust = 0.5))
precip_plot_2022

Soil_plot <- ggplot()+
  geom_line(data = subset(SM_dat, Datetime > as.POSIXct("2023-07-20 00:00") & Datetime < as.POSIXct("2023-10-11 00:00")),aes(x= Datetime, y=soil_water_content),color="cyan3")+
  geom_point(data = daily_predawn, aes(x=Datetime, y=-1*(psi_pd/10)),color="black")+
  geom_line(data = daily_predawn, aes(x=Datetime, y=-1*(psi_pd/10)),color="black")+
  geom_vline(xintercept=as.Date("2023-06-27 00:00:00"),color="grey40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-07-26"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-08-28"),color="gray40",linetype="dashed")+
  scale_y_continuous(
    # Features of the first axis
    name = "30 cm Soil water content (m3/m3)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*10, name="Est. soil water potential (KPa)")
  ) +
  theme_bw()

Soil_plot

ggsave(filename = "LC_2023/2023_weather/weather_graphs/soil_plot.png", plot = Soil_plot, width = 5,height = 2,units = "in",dpi = 300)

#summary of 2022

temp_precip_2022 <- ggplot(subset(weather_daily_2022, Date < "2022-10-15" & Date > "2022-04-11"), aes(x=Date))+
  geom_vline(xintercept=as.Date("2022-06-10"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-07-15"),color="grey40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-08-15"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-09-15"),color="gray40",linetype="dashed")+
  geom_line(aes(y=max_temp, color = "red4"))+
  geom_line(aes(y=daily_precip/5),color="blue3")+
  ggtitle(label = "2022 growing season")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  # Custom the Y scales:
  scale_y_continuous(
    # Features of the first axis
    name = "Daily max T (˚C)",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*5, name="Daily precip (mm)")
  )

temp_precip_2022

##second go...
growth_dat <- read.csv(file= "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")
str(growth_dat)
growth_dat_2022 <- subset(growth_dat, select = c("ID","V299","V335","V357","V385","V421","V497")) %>% summarise(
  V299 = mean(V299),
  V335 = mean(V335),
  V357 = mean(V357),
  V385 = mean(V385),
  V421 = mean(V421),
  V497 = mean(V497)
)
library(reshape2)
library(tidyr)

growth_dat_2022 <- melt(growth_dat_2022, value.name="Volume")
growth_dat_2022$Date <- ymd(c("2022-05-08","2022-06-13","2022-07-05","2022-08-22","2022-09-07","2022-10-15"))

de<-data.frame("V200","2022-04-11","123.0634")
names(de)<-c("variable","Date","Volume")

growth_dat_2022 <- rbind(growth_dat_2022, de)
growth_dat_2022$Volume <- as.numeric(growth_dat_2022$Volume)
growth_dat_2022 <- growth_dat_2022[c(7,1:6),]

#fit a linear fit to growth...
library(lubridate)
ApproxFun2022 <- approxfun(x = growth_dat_2022$Date, y = growth_dat_2022$Volume)
Dates2022 <- seq.Date(ymd("2022-04-12"), ymd("2022-10-14"), by = 1)
LinearFit2022 <- ApproxFun2022(Dates2022)
head(LinearFit2022)

gro_season_2022 <- subset(weather_daily_2022, Date < "2022-10-15" & Date > "2022-04-11")

gro_season_2022$growth <- LinearFit2022

weather_2022_long <- gather(gro_season_2022,measurement,value,max_temp,daily_precip,growth)



weather_2022 <- ggplot((subset(weather_2022_long, Date < "2022-10-15" & Date > "2022-04-11")), aes(x=Date,y=value))+
  geom_line() + theme_bw() + facet_grid(measurement ~.,scales = "free_y",
                                        switch = "y", 
                                        labeller = as_labeller(c(daily_precip = "Daily precip (mm)", growth = "growth (cm3)", max_temp = "Daily max T (˚C)")))+
  ylab(NULL)+
  theme(strip.background = element_blank(), # remove the background
       strip.placement = "outside")+
  geom_vline(xintercept=as.Date("2022-06-10"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-07-15"),color="grey40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-08-15"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2022-09-15"),color="gray40",linetype="dashed")

weather_2022

ggsave(weather_2022,filename="LC_2023/2023_weather/weather_graphs/SPAD_timing_weather.png",height= 5.5, width = 5.5, units = "in",dpi=300)

#2023

str(growth_dat)
growth_dat_2023 <- subset(growth_dat, select = c("ID","V664","V801")) %>% summarise(
  V664 = mean(V664),
  V801 = mean(V801),
)



growth_dat_2023 <- melt(growth_dat_2023, value.name="Volume")
growth_dat_2023$Date <- ymd(c("2023-06-01","2023-10-11"))

de_2023<-data.frame("V500","2022-04-20","1535.833")
names(de_2023)<-c("variable","Date","Volume")

growth_dat_2023 <- rbind(growth_dat_2023, de_2023)
growth_dat_2023$Volume <- as.numeric(growth_dat_2023$Volume)
growth_dat_2023 <- growth_dat_2023[c(3,1,2),]

#fit a linear fit to growth...
library(lubridate)
ApproxFun2023 <- approxfun(x = growth_dat_2023$Date, y = growth_dat_2023$Volume)
Dates2023 <- seq.Date(ymd("2023-04-20"), ymd("2023-10-11"), by = 1)
LinearFit2023 <- ApproxFun2023(Dates2023)
head(LinearFit2023)


gro_season_2023 <- subset(weather_daily, Date < "2023-10-11" & Date > "2023-04-11")

gro_season_2023$growth <- LinearFit2023

weather_2023_long <- gather(gro_season_2023,measurement,value,max_temp,daily_precip,growth,max_VPD)



weather_2023 <- ggplot((subset(weather_2023_long, Date < "2023-10-11" & Date > "2023-04-11")), aes(x=Date,y=value))+
  geom_line() + theme_bw() + facet_grid(measurement ~.,scales = "free_y",
                                        switch = "y", 
                                        labeller = as_labeller(c(daily_precip = "Daily precip (mm)", growth = "growth (cm3)", max_temp = "Daily max T (˚C)")))+
  ylab(NULL)+
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside")+
  geom_vline(xintercept=as.Date("2023-07-15"),color="grey40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-08-15"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-09-15"),color="gray40",linetype="dashed")

weather_2023

ggsave(weather_2023,filename="LC_2023/2023_weather/weather_graphs/SPAD_timing_weather_2023.png",height= 5.5, width = 5.5, units = "in",dpi=300)



weather_2023_long_2 <- gather(gro_season_2023,measurement,value,max_temp,daily_precip,max_VPD)

weather_2023_2 <- ggplot((subset(weather_2023_long_2, Date < "2023-10-11" & Date > "2023-04-11")), aes(x=Date,y=value))+
  geom_line() + theme_bw() + facet_grid(measurement ~.,scales = "free_y",
                                        switch = "y", 
                                        labeller = as_labeller(c(daily_precip = "Daily precip (mm)", max_VPD = "Daily max VPD (kPa)", max_temp = "Daily max T (˚C)")))+
  ylab(NULL)+
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside")+
  geom_vline(xintercept=as.Date("2023-06-27"),color="grey40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-07-26"),color="gray40",linetype="dashed")+
  geom_vline(xintercept=as.Date("2023-08-28"),color="gray40",linetype="dashed")

weather_2023_2             

ggsave(filename = "LC_2023/2023_weather/weather_graphs/weather_diurnal_campaigns.png",plot = weather_2023_2, height = 8, width = 8, units = "in", dpi = 300)
