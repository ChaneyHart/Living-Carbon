#graph 

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)

#thermal data
UAV_thermal <- read.csv("LC_2023/2023_growth_inventory_analysis/processed/UAV_thermal.csv")

UAV_thermal_clean <- UAV_thermal[,-c(1,2,11,13,15,17,19)]

colnames(UAV_thermal_clean) <- c("Row","Column","ID","Construct","Event","area","geometry","temp_0800","temp_0900","temp_1000","temp_1100","temp_1200")

UAV_long <- gather(UAV_thermal_clean,key="Time",value = "temp",temp_0800,temp_0900,temp_1000,temp_1100,temp_1200)
UAV_long$Time <- substr(UAV_long$Time,6,9)
UAV_long$Time <- as.numeric(unlist(UAV_long$Time))
UAV_long$Time <- UAV_long$Time/100
UAV_long$Time <- as.factor(UAV_long$Time)

UAV_long <- UAV_long %>% mutate(Class = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "intermediate",
  Event == "13-15E" | Event == "2H" ~ "high",
  Event == "16-20" | Event == "8-9D" ~ "Control",
  Event == "3" ~ "WT",
  Event == "Wild type"~"Border",
  Event == "1"|Event =="1C"|Event =="13-15B"|Event =="4B"|Event=="5"|Event=="7"~"unanalyzed")
  )


UAV_thermal_summary <- UAV_long %>% group_by(Class,Time) %>% summarize(
  n = n(),
  temp_sd = sd(temp),
  temp = mean(temp),
  temp_se = temp_sd/(sqrt(n))
)

my_colors3 <- c("burlywood1","burlywood1","chartreuse4","gray0","burlywood1","chartreuse4","burlywood1","indianred3","burlywood1","burlywood1","indianred3","indianred3","burlywood1","gray0","burlywood1")
my_colors4 <- c("grey50","gray0","chartreuse4","indianred3","wheat3","grey")

UAV_temp_plot <- ggplot(UAV_thermal_summary, aes(x=Class,y=temp,color=Class))+
  geom_point()+
  geom_errorbar(aes(ymin = temp-temp_se, ymax= temp+temp_se))+
  scale_color_manual(values = my_colors4)+
  facet_grid(.~Time)+
  ylab("Tleaf_UAV (˚C)")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

UAV_temp_plot

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/figures/UAV_temp_plot.png",dpi = 300, width = 8, height = 4, units = "in")

ggplot(subset(UAV_thermal_summary,Time==8), aes(x=Event,y=temp,color=Event))+
  geom_point()+
  geom_errorbar(aes(ymin =temp-temp_se, ymax= temp+temp_se),width=0.2)+
  scale_color_manual(values = my_colors3)+
  theme_bw()


