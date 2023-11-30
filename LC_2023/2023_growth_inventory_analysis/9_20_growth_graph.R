#9_20_23 graph mensuration data

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)



#bring in growth data w/ info about each tree
# Data has been cleaned in 2023_field_season_data_cleaning.R
# Negative growth and missing data was imputted for diameter and height. No outliers were removed
# CT1 was excluded as were trees that were damaged and replaced or died.

#read in data
growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
#check data types
growth_check <- growth %>% group_by(event_short)
group_data(growth_check)
str(growth)
growth$D385 <- as.numeric(growth$D385)

#calculate volume


growth$V49 = ((3.14159*(((growth$D49/10)/2)^2))*(growth$H49/10))
growth$V144 = ((3.14159*(((growth$D144/10)/2)^2))*(growth$H144/10))
growth$V299 = ((3.14159*(((growth$D299/10)/2)^2))*(growth$H299/10))
growth$V335 = ((3.14159*(((growth$D335/10)/2)^2))*(growth$H335/10))
growth$V357 = ((3.14159*(((growth$D357/10)/2)^2))*(growth$H357/10))
growth$V385 = ((3.14159*(((growth$D385/10)/2)^2))*(growth$H385/10))
growth$V421 = ((3.14159*(((growth$D421/10)/2)^2))*(growth$H421/10))
growth$V497 = ((3.14159*(((growth$D497/10)/2)^2))*(growth$H497/10))
growth$V664 = ((3.14159*(((growth$D664/10)/2)^2))*(growth$H664/10))
growth$V801 = ((3.14159*(((growth$D801/10)/2)^2))*(growth$H801/10))

#Height to diameter ratio
growth$HD801 <- growth$H801/growth$D801
growth$HD49 <- growth$H49/growth$V49



#Convert dataframe to long format to look at changes over time for each tree
Height_long <- pivot_longer(growth, cols = c(H49,H144,H299,H335,H357,H385,H421,H497,H664,H801), names_to = "Days", values_to = "Height")
Height_long <- Height_long %>% mutate(Days = as.numeric(substr(Height_long$Days,2,4)))
Height_long$meters <- Height_long$Height/1000

Diam_long <- pivot_longer(growth, cols = c(D49,D144,D299,D335,D357,D385,D421,D497,D664,D801), names_to = "Days", values_to = "Diameter")
Diam_long <- Diam_long %>% mutate(Days = as.numeric(substr(Diam_long$Days,2,4)))

Volume_long <- pivot_longer(growth, cols = c(V49,V144,V299,V335,V357,V385,V421,V497,V664,V801), names_to = "Days", values_to = "Volume")
Volume_long <- Volume_long %>% mutate(Days = as.numeric(substr(Volume_long$Days,2,4)))


#Create a custom color scale
library(RColorBrewer)
myColors <- c('cyan4',
              'red4',
              'gray')


names(myColors) <- levels(growth$construct)
colorScale <- scale_fill_manual(name = "construct",values = myColors)

ht_full <- ggplot(Height_long, aes(x = Days, y = meters))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Height (m)")


ht_full
ggsave("LC_2023/2023_growth_inventory_analysis/ht_timeline_2023.png", plot = ht_full, width = 6, height = 3, units = "in", dpi = 300)

d_full <- ggplot(Diam_long, aes(x = Days, y = Diameter))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("days since planting")+
  ylab("Diameter (mm)")


d_full
ggsave("LC_2023/2023_growth_inventory_analysis/diam_timeline_2023.png", plot = d_full, width = 6, height = 3, units = "in", dpi = 300)


vol_full <- ggplot(Volume_long, aes(x = Days, y = (Volume)))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("days since planting")+
  ylab("Volume (cubic cm)")



vol_full
ggsave("LC_2023/2023_growth_inventory_analysis/vol_timeline_2023.png", plot = vol_full, width = 6, height = 3, units = "in", dpi = 300)

### same thing but broken down by event means######


#Height###

Height_long_II <- Height_long %>% group_by(event_short,Days,block)

Height_long_II_summary <- Height_long_II %>% dplyr::summarise(
  height = mean(meters),
  n = n(),
  height_sd = sd(meters),
  height_se = height_sd/(sqrt(n))
)

Height_long_III <- Height_long %>% group_by(event_short,Days)

Height_long_III_summary <- Height_long_III %>% dplyr::summarise(
  height = mean(meters),
  n = n(),
  height_sd = sd(meters),
  height_se = height_sd/(sqrt(n))
)

ht_full_event <- ggplot(Height_long_III_summary, aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Height (m)")
#large block
ht_full_event_lb <- ggplot(subset(Height_long_II_summary, block == "large"), aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("height (m)")
##small block
ht_full_event_sb <- ggplot(subset(Height_long_II_summary, block == "small"), aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("height (m)")



ht_full_event
ht_full_event_lb
ht_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis//ht_timeline_event_2023.png", plot = ht_full_event, width = 6, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_lb.png", plot = ht_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_sb.png", plot = ht_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)



###diameter by event #####

diam_long_II <- Diam_long %>% group_by(event_short,Days, block)

diam_long_II_summary <- diam_long_II %>% dplyr::summarise(
  diam = mean(Diameter),
  n = n(),
  diameter_sd = sd(Diameter),
  diameter_se = diameter_sd/(sqrt(n))
)

diam_long_III <- Diam_long %>% group_by(event_short,Days)

diam_long_III_summary <- diam_long_III %>% dplyr::summarise(
  diam = mean(Diameter),
  n = n(),
  diameter_sd = sd(Diameter),
  diameter_se = diameter_sd/(sqrt(n))
)


#all
Diam_full_event <- ggplot(diam_long_III_summary, aes(x = Days, y = diam))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

#large block
Diam_full_event_lb <- ggplot(subset(diam_long_II_summary, block == "large"), aes(x = Days, y = diam))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")
##small block
Diam_full_event_sb <- ggplot(subset(diam_long_II_summary, block == "small"), aes(x = Days, y = diam))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

Diam_full_event
Diam_full_event_lb
Diam_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis//Diam_timeline_event_2023.png", plot = Diam_full_event, width = 6, height = 3, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_lb.png", plot = Diam_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_sb.png", plot = Diam_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)


#####volume by event ###########

volume_long_II <- Volume_long %>% group_by(event_short,Days, block)
volume_long_III <- Volume_long %>% group_by(event_short,Days)

volume_long_II_summary <- volume_long_II %>% dplyr::summarise(
  vol = mean(Volume),
  n = n(),
  vol_sd = sd(Volume),
  vol_se = vol_sd/(sqrt(n))
)

volume_long_III_summary <- volume_long_III %>% dplyr::summarise(
  vol = mean(Volume),
  n = n(),
  vol_sd = sd(Volume),
  vol_se = vol_sd/(sqrt(n))
)


vol_full_event <- ggplot(volume_long_III_summary, aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")+
  theme_bw()

vol_full_event
#large block
vol_full_event_lb <- ggplot(subset(volume_long_II_summary, block == "large"), aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")

#small block

vol_full_event_sb <- ggplot(subset(volume_long_II_summary, block == "small"), aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")


vol_full_event
vol_full_event_lb
vol_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis/vol_timeline_event_2023.png", plot = vol_full_event, width = 6, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_lb.png", plot = vol_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_sb.png", plot = vol_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)


####table to summarize and graph data
#mean volume per event

Volume_summary <- data_frame()
event_volume <- as.data.frame(aggregate(V801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(V801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(V801 ~ event_short, data = growth, FUN = length)

event_volume_summary <- inner_join(event_volume, event_var, by = "event_short")
event_volume_summary <- inner_join(event_volume_summary, event_n, by = "event_short")

colnames(event_volume_summary) <- c("Event", "mean_Volume","var","n")

event_volume_summary$se <- sqrt(event_volume_summary$var/event_volume_summary$n)
event_volume_summary$tst <- qt(0.975, event_volume_summary$n)
event_volume_summary$sd <- sqrt(event_volume_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_volume_summary$upper_volume <- event_volume_summary$mean_Volume + (event_volume_summary$tst*event_volume_summary$se)
event_volume_summary$lower_volume <- event_volume_summary$mean_Volume - (event_volume_summary$tst*event_volume_summary$se)


write.csv(event_volume_summary, file = "LC_2023/2023_growth_inventory_analysis/event_volume_table.csv")
event_volume_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

event_volume_summary <- event_volume_summary %>% mutate(tier = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "top",
  Event == "13-15E" | Event == "2H" ~ "poor",
  Event == "16-20" | Event == "8-9D" ~ "Control",
  Event == "CT3" ~ "WT",
  Event == "13-15B" | Event == "1" | Event == "1C" | Event == "7"| Event == "4B"| Event == "5" ~ "unanalyzed"))


myColors4 <- c('gray4',
              'chartreuse4',
              'indianred3',
              'blue3',
              'gray')


names(myColors4) <- levels(event_volume_summary$tier)
colorScale <- scale_color_manual(name = "tier",values = myColors)


volume_plot <- ggplot(event_volume_summary, aes(x=reorder(Event,mean_Volume),y=mean_Volume, color = tier))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_volume,ymax=upper_volume))+
  xlab("Event")+
  ylab("Volume index (cm^3)")+
  scale_color_manual(values = myColors4)+
  theme_bw()

#volume_plot_2 <- volume_plot + colorScale
volume_plot

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/volume_plot_2023.png", plot = volume_plot, width = 8, height =5, units= "in",dpi = 300)

###Diameter summary


diameter_summary <- data_frame()
event_diameter <- as.data.frame(aggregate(D801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(D801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(D801 ~ event_short, data = growth, FUN = length)

event_diameter_summary <- inner_join(event_diameter, event_var, by = "event_short")
event_diameter_summary <- inner_join(event_diameter_summary, event_n, by = "event_short")

event_diameter_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")
event_diameter_summary <- event_diameter_summary %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))


myColors <- c('gray',
              'gray0',
              'indianred3',
              'chartreuse4')


names(myColors) <- levels(event_diameter_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)

colnames(event_diameter_summary) <- c("Event", "mean_diameter","var","n","construct","tier")

event_diameter_summary$se <- sqrt(event_diameter_summary$var/event_diameter_summary$n)
event_diameter_summary$tst <- qt(0.975, event_diameter_summary$n)
event_diameter_summary$sd <- sqrt(event_diameter_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_diameter_summary$upper_diameter <- event_diameter_summary$mean_diameter + (event_diameter_summary$tst*event_diameter_summary$se)
event_diameter_summary$lower_diameter <- event_diameter_summary$mean_diameter - (event_diameter_summary$tst*event_diameter_summary$se)

write.csv(event_diameter_summary, file = "LC_2023/2023_growth_inventory_analysis/event_diameter_table.csv")

d_plot<- ggplot(event_diameter_summary, aes(x=reorder(Event,mean_diameter),y=mean_diameter,color= tier))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_diameter,ymax=upper_diameter))+
  xlab("Event")+
  ylab("diameter (mm)")+
  scale_color_manual(values = myColors4)+
  theme_bw()

d_plot 

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/final_diam_plot.png",plot=d_plot, width = 8, height = 5, units = "in",dpi=300)
##Height summary#####
height_summary <- data_frame()
event_height <- as.data.frame(aggregate(H801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(H801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(H801 ~ event_short, data = growth, FUN = length)

event_height_summary <- inner_join(event_height, event_var, by = "event_short")
event_height_summary <- inner_join(event_height_summary, event_n, by = "event_short")

colnames(event_height_summary) <- c("Event", "mean_height","var","n")

event_height_summary$se <- sqrt(event_height_summary$var/event_height_summary$n)
event_height_summary$tst <- qt(0.975, event_height_summary$n)
event_height_summary$sd <- sqrt(event_height_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_height_summary$upper_height <- event_height_summary$mean_height + (event_height_summary$tst*event_height_summary$se)
event_height_summary$lower_height <- event_height_summary$mean_height - (event_height_summary$tst*event_height_summary$se)

write.csv(event_height_summary, file = "LC_2023/2023_growth_inventory_analysis/event_height_table.csv")

event_height_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

event_height_summary <- event_height_summary %>% mutate(tier = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "top",
  Event == "13-15E" | Event == "2H" ~ "poor",
  Event == "16-20" | Event == "8-9D" ~ "Control",
  Event == "CT3" ~ "WT",
  Event == "13-15B" | Event == "1" | Event == "1C" | Event == "7"| Event == "4B"| Event == "5" ~ "unanalyzed"))


names(myColors) <- levels(event_height_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)

ht_plot <- ggplot(event_height_summary, aes(x=reorder(Event,mean_height),y=mean_height, color=tier))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_height,ymax=upper_height))+
  xlab("Event")+
  ylab("height (mm)")+
  scale_color_manual(values = myColors4)+
  theme_bw()


ht_plot

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/final_ht_plot.png",plot= ht_plot, width= 8, height = 5, units = "in",dpi=300)

####
#Height to diameter ratio

HD_summary <- data_frame()
event_HD <- as.data.frame(aggregate(HD801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(HD801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(HD801 ~ event_short, data = growth, FUN = length)

event_HD_summary <- inner_join(event_HD, event_var, by = "event_short")
event_HD_summary <- inner_join(event_HD_summary, event_n, by = "event_short")

colnames(event_HD_summary) <- c("Event", "mean_HD","var","n")

event_HD_summary$se <- sqrt(event_HD_summary$var/event_HD_summary$n)
event_HD_summary$tst <- qt(0.975, event_HD_summary$n)

event_HD_summary$sd <- sqrt(event_HD_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_HD_summary$upper_HD <- event_HD_summary$mean_HD + (event_HD_summary$tst*event_HD_summary$se)
event_HD_summary$lower_HD <- event_HD_summary$mean_HD - (event_HD_summary$tst*event_HD_summary$se)

write.csv(event_HD_summary, file = "LC_2023/2023_growth_inventory_analysis/event_HD_table.csv")

event_HD_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

event_HD_summary <- event_HD_summary %>% mutate(tier = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "top",
  Event == "13-15E" | Event == "2H" ~ "poor",
  Event == "16-20" | Event == "8-9D" ~ "Control",
  Event == "CT3" ~ "WT",
  Event == "13-15B" | Event == "1" | Event == "1C" | Event == "7"| Event == "4B"| Event == "5" ~ "unanalyzed"))


names(myColors) <- levels(event_HD_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)

HD_plot <- ggplot(event_HD_summary, aes(x=reorder(Event,mean_HD),y=mean_HD, color=tier))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_HD,ymax=upper_HD))+
  xlab("Event")+
  ylab("hieght to diameter ratio (mm)")+
  scale_color_manual(values = myColors4)+
  theme_bw()


HD_plot

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/final_HD_plot.png",plot= HD_plot, width= 8, height = 5, units = "in",dpi=300)


#Explore spatial variation in covariates
#read in covariates
covariates_inner <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/covariates_inner.csv")
covariates_full <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/covariates_full.csv")


covariates_set_inner <- inner_join(growth, covariates_inner, by = "ID")
covariates_set_full <- full_join(growth, covariates_full, by="ID")

####SPAD

str(covariates_set_full)

ggplot(subset(covariates_set_full, avg_SPAD < 100), aes(x=event, y=avg_SPAD))+
  geom_point()

#spatially
library(viridisLite)
library(viridis)
library(latticeExtra)

levelplot(avg_SPAD ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))



levelplot(SPAD_June_2022 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_July_2022 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_aug_2022 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_sep_2022 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_July_2023 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_Aug_2023 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(SPAD_Sep_2023 ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))



##NDVI and other indices

str(covariates_full)
levelplot(GNDVI_mean ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

ggplot(subset(covariates_set_full, avg_SPAD < 100), aes(x=GNDVI_mean, y=SPAD_aug_2022))+
  geom_point()

levelplot(GRVI_mean ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(GCI_mean ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(ARI_mean ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))

levelplot(RECI_mean ~-column*row, subset(covariates_set_full,avg_SPAD < 100),panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))


#crown area via shadows
#scatter plot

#log transformed
ggplot(UAV_join, aes(x=log(V497),y=log(shadow_area), color = HD))+
  geom_point()+
  scale_color_viridis()+
  geom_abline(intercept=-4.2078,slope = 0.4972)


#would expect that trees with a lower height to diameter ratio may be branchier and tend to be above the fitted line.
?geom_abline

branchiness_lm <- lm(log(shadow_area)~log(V497), data =UAV_join )
summary(branchiness_lm)
plot(branchiness_lm, which = 1)
#outliers to look at to see if they are branchier than expected

#residuals from this relationship could feasibly be related to deviations of crown area from stem volume due to branchiness

ggplot(UAV_join, aes(y=HD, x=branchiness_lm$residuals))+
  geom_point()

#there may be a slight relationship here but HD ratio does not seem to explain branchiness residuals

#how to branchiness residuals vary spatially?

library(latticeExtra)

levelplot(branchiness_lm$residuals ~-column*row, UAV_join,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
#this may be helpful. It does seem to match with where I see branchier trees.

levelplot(shadow_area ~-column*row, UAV_join,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))


