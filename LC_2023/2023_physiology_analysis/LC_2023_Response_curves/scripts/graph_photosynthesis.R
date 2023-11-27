##Plot A-Ci graphs
library(ggplot2)
library(dplyr)
library(lubridate)

#import cleaned ACI data

ACI_summary_tree <- read.csv("LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_summary_tree_cleaned.csv")
ACI_summary_tree$PhiPS2 <- as.numeric(ACI_summary_tree$PhiPS2)
ACI_summary_tree$ETR <- as.numeric(ACI_summary_tree$ETR)

#import photosynthetic parameters
ACI_parameters <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/aci_parameters_list.csv")

#import event summary

ACI_event_summary <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_event_summary.csv")

#light_response_curves

#import weather data for graphing
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

#tier summary
ACI_tier_summary <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_tier_summary.csv")
ACI_tier_summary_date <- read.csv(file ="LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_tier_summary_date.csv")

####graph_sample dates with temp and VPD

Sample_date_plot <- ggplot()+
  geom_line(data = subset(weather_daily, Datetime > as.POSIXct("2023-06-01 00:00") & Datetime < as.POSIXct("2023-9-01 00:00")),aes(x= Datetime, y=max_temp),color="red3")+
  geom_line(data = subset(weather_daily, Datetime > as.POSIXct("2023-06-01 00:00") & Datetime < as.POSIXct("2023-9-01 00:00")),aes(x= Datetime, y=max_VPD*10),color="blue3")+
  geom_point(aes(x=as.POSIXct("2023-06-14 12:00"),y=0))+
  geom_point(aes(x=as.POSIXct("2023-06-22 12:00"),y=0))+
  geom_point(aes(x=as.POSIXct("2023-07-18 12:00"),y=0))+
  geom_point(aes(x=as.POSIXct("2023-07-22 12:00"),y=0))+
  geom_point(aes(x=as.POSIXct("2023-07-29 12:00"),y=0))+
  geom_point(aes(x=as.POSIXct("2023-08-29 12:00"),y=0))+
  scale_y_continuous(
    # Features of the first axis
    name = "max air temp (˚C)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./10, name="max VPD (KPa)")
  ) +
  theme_bw()

Sample_date_plot
ggsave(filename = "LC_2023/2023_weather/weather_graphs/Sample_date_weather.png",plot = Sample_date_plot, dpi = 300)
#### graph summary ACI curves

library(RColorBrewer)
my_colors <- c("gray0","chartreuse4","indianred3","gray")
my_colors2 <- c("chartreuse4","gray0","chartreuse4","indianred3","indianred3","indianred3","gray0","gray0")


ACI_tier_plot <- ggplot(ACI_tier_summary, aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_line(linetype=2,alpha=0.5)+
  geom_errorbar(aes(ymax=A_upper,ymin=A_lower))+
  scale_color_manual(values = my_colors)+
  theme_bw()

ACI_tier_plot
ggsave(filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/A_Ci_tier_summary.png",plot=ACI_tier_plot, dpi=300)

ggplot(ACI_event_summary, aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(linewidth = 0.5,se=TRUE )+
  geom_errorbar(aes(ymax=A_upper,ymin=A_lower))+
  scale_color_manual(values= my_colors2)


#Graph A-CI curves dependent on sample date


ggplot(subset(ACI_tier_summary_date, sample_date.x == "Jun14"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)

ggplot(subset(ACI_tier_summary_date, sample_date.x == "Jun22"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)

ggplot(subset(ACI_tier_summary_date, sample_date.x == "July18"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)
#2H missing here

ggplot(subset(ACI_tier_summary_date, sample_date.x == "July20"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)

#2H and 5A missing here
ggplot(subset(ACI_tier_summary_date, sample_date.x == "July27"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)

#2H missing here
ggplot(subset(ACI_tier_summary_date, sample_date.x == "Aug29"), aes(x=Ci,y=A, color = tier))+
  geom_point()+
  geom_smooth(se=FALSE,linewidth=0.7)+
  geom_errorbar(aes(ymin=A_lower,ymax=A_upper))+
  scale_color_manual(values = my_colors)






