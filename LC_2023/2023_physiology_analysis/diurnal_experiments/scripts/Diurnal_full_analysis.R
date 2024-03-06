##Combined analysis of diurnal experiments

library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(tidyr)
library(stringi)
library(janitor)

#June

June_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_tree.csv)")
June_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_tree.csv)")

June_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_event.csv)")
June_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_event.csv)")

June_summary_event <- inner_join(June_summary_6800_event, June_summary_600_event, by = c("event_short","Timepoint"))

June_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_class.csv)")
June_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_class.csv)")

June_summary_class <- inner_join(June_summary_6800_class, June_summary_600_class, by = c("Class","Timepoint"))
June_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_28_fluor_set_diurnal_summary_class.csv))")


#July

July_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_tree.csv)")
July_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_tree.csv)")

July_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_event.csv)")
July_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_event.csv)")

July_summary_event <- inner_join(July_summary_6800_event, July_summary_600_event, by = c("event_short","Timepoint","treatment"))

July_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_class.csv)")
July_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_class.csv)")

July_summary_class <- inner_join(July_summary_6800_class, July_summary_600_class, by = c("Class","Timepoint","treatment"))
July_fluor_set_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_event.csv)")
July_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_class.csv)")

#August

August_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_tree.csv)")
August_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_tree.csv)")
August_predawn_summary_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_tree.csv)")

August_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_event.csv)")
August_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_event.csv)")
August_predawn_summary_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_event.csv)")

August_summary_event <- inner_join(August_summary_6800_event, August_summary_600_event, by = c("event_short","Timepoint","treatment"))

August_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_class.csv)")
August_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_class.csv)")
August_predawn_summary_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_class.csv)")

August_summary_class <- inner_join(August_summary_6800_class, August_summary_600_class, by = c("Class","Timepoint","treatment"))

August_fluor_set_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_event.csv)")
August_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_class.csv)")


#######graph######


#define colors
library(RColorBrewer)
my_colors <- c("gray0","chartreuse4","indianred3","grey")

library(ggsignif)
library(gridExtra)
June_assimilation <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  xlab("Hour")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)+
  theme_bw()


June_assimilation

June_gsw <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()
June_gsw

June_ETR <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()
June_ETR

June_PhiPS2 <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PhiPS2), size = 2)+
  geom_line(aes(y = PhiPS2))+ 
  geom_errorbar(aes(ymin = PhiPS2-PhiPS2_se, ymax = PhiPS2+PhiPS2_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("efficiency of PSII")+
  xlab("Hour")+
  theme_bw()

June_PhiPS2

June_Tleaf <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = Tleaf), size = 2)+
  geom_line(aes(y = Tleaf))+ 
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")+
  xlab("Hour")+
  theme_bw()
June_Tleaf

June_diurnal_phys_plot <- grid.arrange(June_assimilation, June_ETR, June_gsw,June_PhiPS2,June_Tleaf,nrow=2)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_diurnal_phys_plot.png", plot = June_diurnal_phys_plot, height = 8, width = 12, units = "in", dpi = 300)

#photorespiration phys

png(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/A_ETR.png",height = 4, width = 4, units = "in",res = 300)
ggplot(June_fluor_set, aes(x=A-R,y=ETR,color = Class))+
  geom_point()+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Gross assimilation (µmol/m^2*s)")+
  theme_bw()

dev.off()


June_Jo <- ggplot(June_fluor_set_summary_class, aes(x = Timepoint, group=Class,color = Class))+
  geom_point(aes(y = Jo), size = 2)+
  geom_line(aes(y = Jo))+ 
  geom_errorbar(aes(ymin = Jo-Jo_se, ymax = Jo+Jo_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")
June_Jo

June_Rl <- ggplot(June_fluor_set_summary_class, aes(x = Timepoint, group=Class,color = Class))+
  geom_point(aes(y = Rl), size = 2)+
  geom_line(aes(y = Rl))+ 
  geom_errorbar(aes(ymin = Rl-Rl_se, ymax = Rl+Rl_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Rl (µmol CO2 /m^2*s)")
June_Rl


June_Jo_ratio <- ggplot(subset(June_fluor_set_summary_class, Timepoint != 18), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  geom_line(aes(y = Jo_ratio))+ 
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio+Jo_ratio_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
June_Jo_ratio


June_photorespiration_plot <- grid.arrange(June_Jo, June_Rl, June_Jo_ratio, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_photorespiration_plot.png",plot = June_photorespiration_plot, width = 8, height = 6, dpi = 300)


library(RColorBrewer)
my_colors <- c("gray0","chartreuse4","indianred3","grey")
my_colors2 <- c("chartreuse4","gray0","chartreuse4","indianred3","indianred3","indianred3","gray0","gray")

library(ggsignif)
library(gridExtra)
July_assimilation <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  ylab("Assimilation (µmol/m^2*s)")+
  scale_color_manual(values = my_colors)+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)+
  theme_bw()

July_assimilation

July_gsw <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")+
  theme_bw()
July_gsw

July_summary_class$WUE <- July_summary_class$A/July_summary_class$gsw

July_WUE <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = WUE), size = 2)+
  geom_line(aes(y = WUE, linetype = treatment))+
  scale_color_manual(values = my_colors)+
  ylab("Water use efficiency")
July_WUE


July_ETR <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  theme_bw()
July_ETR

July_diurnal_phys_plot <- grid.arrange(July_assimilation, July_gsw, July_ETR, nrow=3)

July_Jo <- ggplot(subset(July_fluor_set_class), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo), size = 2)+
  geom_line(aes(y = Jo, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo-Jo_se, ymax = Jo+Jo_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")
July_Jo

July_Jo_ratio <- ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  geom_line(aes(y = Jo_ratio, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio+Jo_ratio_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
July_Jo_ratio

#gsw curves

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation")

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("ETR")




#August
str(August_summary_class)
August_summary_class$Timepoint <- as.numeric(August_summary_class$Timepoint)

August_assimilation_class <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A,linetype = treatment,linetype = treatment))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)

August_assimilation_class

August_assimilation_event <- ggplot(August_summary_event, aes(x = Timepoint, color = event_short, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A,linetype = treatment,linetype = treatment))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors2)+
  ylab("Assimilation (µmol/m^2*s)")

August_assimilation_event

August_gsw <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw,linetype = treatment))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")
August_gsw

August_ETR <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")
August_ETR

August_diurnal_phys_plot <- grid.arrange(August_assimilation_class, August_gsw, August_ETR, nrow=3)

August_Jo <- ggplot(August_fluor_set_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo), size = 2)+
  geom_line(aes(y = Jo,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo-Jo_se, ymax = Jo+Jo_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")
August_Jo

August_Jo_ratio <- ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  geom_line(aes(y = Jo_ratio,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio+Jo_ratio_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
August_Jo_ratio

#gsw curves

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation")

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("ETR")


ggplot(subset(August_fluor_set_class, Class = Elite), aes(x = Timepoint, color = treatment))+
  geom_point(aes(y=A))+
  geom_line(aes(y=A))+
  geom_point(aes(y=ETR))+
  geom_line(aes(y=ETR))

ggplot(subset(August_fluor_set_class, Class = Elite), aes(x = Timepoint, color = treatment))+
  geom_point(aes(y=Jo))+
  geom_line(aes(y=Jo))+
  geom_point(aes(y=Jc))+
  geom_line(aes(y=Jc))

#weather dat

library(lubridate)
weather_6_27 <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
str(weather_6_27)
weather_6_27$Datetime <- ymd_hms(weather_6_27$Datetime)
weather_6_27$Date <- as.Date(weather_6_27$Datetime)
weather_6_27$Date
weather_6_27 <- subset(weather_6_27, Date == "2023-06-27")
weather_6_27$hour <- hour(weather_6_27$Datetime)


June_solar_rad <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = solar_radiation_mean), size = 2)+
  geom_line(aes(y = solar_radiation_mean))+ 
  ylab("mean solar radiation (W/m2)")+
  theme_bw()
June_solar_rad

June_air_temp <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = air_temp_mean), size = 2)+
  geom_line(aes(y = air_temp_mean))+ 
  ylab("mean temp (˚C)")+
  theme_bw()
June_air_temp

June_air_VPD <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = VPD_mean), size = 2)+
  geom_line(aes(y = VPD_mean))+ 
  ylab("VPD (kPa)")+
  theme_bw()
June_air_VPD

June_weather_station <- grid.arrange(June_solar_rad, June_air_temp, June_air_VPD, nrow=3)
ggsave(plot = June_weather_station,filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_weather_station.png",height = 8, width = 8, units = "in",dpi = 300 )

#July
library(lubridate)
weather_7_26 <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
str(weather_7_26)
weather_7_26$Datetime <- ymd_hms(weather_7_26$Datetime)
weather_7_26$Date <- as.Date(weather_7_26$Datetime)
weather_7_26$Date
weather_7_26 <- subset(weather_7_26, Date == "2023-07-26")
weather_7_26$hour <- hour(weather_7_26$Datetime)


July_solar_rad <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = solar_radiation_mean), size = 2)+
  geom_line(aes(y = solar_radiation_mean))+ 
  ylab("mean solar radiation (W/m2)")+
  theme_bw()
July_solar_rad

July_air_temp <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = air_temp_mean), size = 2)+
  geom_line(aes(y = air_temp_mean))+ 
  ylab("mean temp (˚C)")+
  theme_bw()
July_air_temp

July_air_VPD <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = VPD_mean), size = 2)+
  geom_line(aes(y = VPD_mean))+ 
  ylab("VPD (kPa)")+
  theme_bw()
July_air_VPD

July_weather_station <- grid.arrange(July_solar_rad, July_air_temp, July_air_VPD, nrow=3)
ggsave(plot = July_weather_station,filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_weather_station.png",height = 8, width = 8, units = "in",dpi = 300 )


