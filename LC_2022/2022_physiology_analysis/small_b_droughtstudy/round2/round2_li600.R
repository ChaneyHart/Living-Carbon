#drought diurnal curves - round 2

install.packages("lubridate")
install.packages("scales")
install.packages("ggthemes", dependencies = TRUE)
library(ggplot2)
library(lubridate)
library(gridExtra)
library((ggthemes))
install.packages("ggsignif")
library(ggsignif)
library(dplyr)




dat_rd2 <- read.csv("small_b_droughtstudy/round2/round2_diurnal_curves_proc.csv")
dat_rd2 <- subset(dat_rd2, select = c(Time,time.hour.,ID,gsw,VPDleaf,PhiPS2,ETR,Tleaf,Qamb))
dat_meta <- read.csv("DBH_H_timeline_CT1_excluded_9_22.csv")
dat_meta <- subset(dat_meta, select = c(row,column,ID,event_short, construct, construct2,block,H419))



dat_rd2 <- inner_join(dat_rd2, dat_meta, by = "ID")

#dat_rd1$Event <- as.factor(dat_rd1$Event)

library(RColorBrewer)
myColors <- c('gray',
              'blue',
              'red'
)
myColors2 <- c('gray', 'red')

names(myColors) <- levels(dat_rd2$Construct2)
colorScale <- scale_fill_manual(name = "Construct2",values = myColors)


#10AM

dat_10 <- subset(dat_rd2, time.hour. =="10")

gsw_10_plot <- ggplot(dat_10, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (µmol/mol s)")+
  ggtitle("10:00")+
  theme(legend.position = "None")

rd2_gsw_10_plot2 <- ggplot(dat_10, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("10:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.6)+
  theme(legend.position = "None")

gsw_10_plot <- gsw_10_plot + colorScale
gsw_10_plot
rd2_gsw_10_plot2 
ggsave("small_b_droughtstudy/round2/gsw_10_rd2.png", plot = rd2_gsw_10_plot2,dpi = 300)

Tleaf_10_plot <- ggplot(dat_10, aes(x = Construct2, y = Tleaf, fill = Construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("Tleaf ˚C")+
  ggtitle("10:00")+
  geom_signif(comparisons = list(c("Grape","Control")),map_signif_level=TRUE, y_position = 34)+
  geom_signif(comparisons = list(c("Poplar","Control")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("Poplar","Grape")),map_signif_level=TRUE, y_position = 35)+
  theme(legend.position = "None")


Tleaf_10_plot <- Tleaf_10_plot + colorScale
Tleaf_10_plot
ggsave("Tleaf_10.png", plot = Tleaf_10_plot, dpi = 300)

#12
#need to average values for consecutive measurements on each tree

dat_12 <- subset(dat_rd2, time.hour. =="12")

dat_12_2 <- inner_join(dat_12_2, dat_meta, by = "ID")

gsw_12_plot <- ggplot(dat_12, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("12:00")+
  theme(legend.position = "None")

rd2_gsw_12_plot2 <- ggplot(dat_12, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("12:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.55)+
  theme(legend.position = "None")

gsw_12_plot <- gsw_12_plot + colorScale
rd2_gsw_12_plot2

gsw_12_plot
ggsave("small_b_droughtstudy/round2/gsw_12_rd2.png", plot = rd2_gsw_12_plot2, dpi = 300)

Tleaf_12_plot <- ggplot(dat_12_2, aes(x = Construct2, y = Tleaf, fill = Construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("Tleaf ˚C")+
  ggtitle("12:00")+
  geom_signif(comparisons = list(c("Grape","Control")),map_signif_level=TRUE, y_position = 30)+
  geom_signif(comparisons = list(c("Poplar","Control")),map_signif_level=TRUE, y_position = 34)+
  geom_signif(comparisons = list(c("Poplar","Grape")),map_signif_level=TRUE, y_position = 33)+
  theme(legend.position = "None")


Tleaf_12_plot <- Tleaf_12_plot + colorScale

ggsave("Tleaf_12.png", plot = Tleaf_12_plot, dpi = 300)
# 14

dat_14 <- subset(dat_rd2, time.hour. =="14")


dat_14_2 <- inner_join(dat_14_2, dat_meta, by = "ID")

gsw_14_plot <- ggplot(dat_14, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("14:00")+
  theme(legend.position = "None")

rd2_gsw_14_plot2 <- ggplot(dat_14, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("14:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.35)+
  theme(legend.position = "None")


gsw_14_plot <- gsw_14_plot + colorScale
gsw_14_plot
rd2_gsw_14_plot2
ggsave("small_b_droughtstudy/round2/gsw_14_rd2.png", plot = rd2_gsw_14_plot2, dpi = 300)

Tleaf_14_plot <- ggplot(dat_14_2, aes(x = Construct2, y = Tleaf, fill = Construct))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("Tleaf ˚C")+
  ggtitle("14:00")+
  geom_signif(comparisons = list(c("Grape","Control")),map_signif_level=TRUE, y_position = 34)+
  geom_signif(comparisons = list(c("Poplar","Control")),map_signif_level=TRUE, y_position = 36)+
  geom_signif(comparisons = list(c("Poplar","Grape")),map_signif_level=TRUE, y_position = 35)+
  theme(legend.position = "None")


Tleaf_14_plot <- Tleaf_14_plot + colorScale

ggsave("Tleaf_14.png", plot = Tleaf_14_plot, dpi = 300)

#16

dat_16 <- subset(dat_rd2, time.hour. =="16")


gsw_16_plot <- ggplot(dat_16, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("16:00")+
  theme(legend.position = "None")

rd2_gsw_16_plot2 <- ggplot(dat_16, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("16:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.4)+
  theme(legend.position = "None")

gsw_16_plot <- gsw_16_plot + colorScale

rd2_gsw_16_plot2
gsw_16_plot

ggsave("small_b_droughtstudy/round2/gsw_16_rd2.png", plot = rd2_gsw_16_plot2, dpi = 300)

Tleaf_16_plot <- ggplot(dat_16_2, aes(x = Construct2, y = Tleaf, fill = Construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("Tleaf ˚C")+
  ggtitle("16:00")+
  geom_signif(comparisons = list(c("Grape","Control")),map_signif_level=TRUE, y_position = 34)+
  geom_signif(comparisons = list(c("Poplar","Control")),map_signif_level=TRUE,y_position = 36)+
  geom_signif(comparisons = list(c("Poplar","Grape")),map_signif_level=TRUE, y_position = 35)+
  theme(legend.position = "None")


Tleaf_16_plot <- Tleaf_16_plot + colorScale

ggsave("Tleaf_16.png", plot = Tleaf_16_plot, dpi = 300)
#############
#stats
lm_10_gsw <- lm(gsw ~ construct2, dat_10)
aov_10_gsw <- aov(gsw ~ construct2, dat_10)
lm_10_Tleaf <- lm(Tleaf ~ construct2, dat_10)
aov_10_Tleaf <- aov(Tleaf ~ construct2, dat_10)
anova(lm_10_gsw)
TukeyHSD(aov_10_gsw)
anova(lm_10_Tleaf)
TukeyHSD(aov_10_Tleaf)

###
lm_12_gsw <- lm(gsw ~ construct2, dat_12)
aov_12_gsw <- aov(gsw ~construct2, dat_12)
lm_12_Tleaf <- lm(Tleaf ~ construct2, dat_12)
aov_12_Tleaf <- aov(Tleaf ~construct2, dat_12)
anova(lm_12_gsw)
TukeyHSD(aov_12_gsw)
anova(lm_12_Tleaf)
TukeyHSD(aov_12_Tleaf)

####
lm_14_gsw <- lm(gsw ~ construct2, dat_14)
aov_14_gsw <- aov(gsw ~construct2, dat_14)
lm_14_Tleaf <- lm(Tleaf ~ construct2, dat_14)
aov_14_Tleaf <- aov(Tleaf ~construct2, dat_14)
anova(lm_14_gsw)
TukeyHSD(aov_14_gsw)
anova(lm_14_Tleaf)
TukeyHSD(aov_14_Tleaf)

#######
lm_16_gsw <- lm(gsw ~ construct2, dat_16)
aov_16_gsw <- aov(gsw ~construct2, dat_16)
lm_16_Tleaf <- lm(Tleaf ~ construct2, dat_16)
aov_16_Tleaf <- aov(Tleaf ~construct2, dat_16)
anova(lm_16_gsw)
TukeyHSD(aov_16_gsw)
anova(lm_16_Tleaf)
TukeyHSD(aov_16_Tleaf)

################

#stats
lm_10_gsw <- lm(gsw ~ event, dat_10)
aov_10_gsw <- aov(gsw ~ event, dat_10)
lm_10_Tleaf <- lm(Tleaf ~ event, dat_10)
aov_10_Tleaf <- aov(Tleaf ~ event, dat_10)
anova(lm_10_gsw)
TukeyHSD(aov_10_gsw)
anova(lm_10_Tleaf)
TukeyHSD(aov_10_Tleaf)

###
lm_12_gsw <- lm(gsw ~ event, dat_12)
aov_12_gsw <- aov(gsw ~event, dat_12)
lm_12_Tleaf <- lm(Tleaf ~ event, dat_12)
aov_12_Tleaf <- aov(Tleaf ~event, dat_12)
anova(lm_12_gsw)
TukeyHSD(aov_12_gsw)
anova(lm_12_Tleaf)
TukeyHSD(aov_12_Tleaf)

####
lm_14_gsw <- lm(gsw ~ event, dat_14)
aov_14_gsw <- aov(gsw ~event, dat_14)
lm_14_Tleaf <- lm(Tleaf ~ event, dat_14)
aov_14_Tleaf <- aov(Tleaf ~event, dat_14)
anova(lm_14_gsw)
TukeyHSD(aov_14_gsw)
anova(lm_14_Tleaf)
TukeyHSD(aov_14_Tleaf)

#######
lm_16_gsw <- lm(gsw ~ event, dat_16)
aov_16_gsw <- aov(gsw ~event, dat_16)
lm_16_Tleaf <- lm(Tleaf ~ event, dat_16)
aov_16_Tleaf <- aov(Tleaf ~event, dat_16)
anova(lm_16_gsw)
TukeyHSD(aov_16_gsw)
anova(lm_16_Tleaf)
TukeyHSD(aov_16_Tleaf)


gsw_plot <- grid.arrange(gsw_10_plot,gsw_12_plot, gsw_14_plot,gsw_16_plot, nrow=2)
ggsave("gsw_plot.png", plot = gsw_plot,height = 5.0,width = 10.0, dpi = 300)

rd2_gsw_plot2 <- grid.arrange(rd2_gsw_10_plot2,rd2_gsw_12_plot2, rd2_gsw_14_plot2,rd2_gsw_16_plot2, nrow=2)
ggsave("small_b_droughtstudy/round2/rd_2_gsw_plot2.png", plot = rd2_gsw_plot2,height = 5.0, dpi = 300)

Tleaf_plot <- grid.arrange(Tleaf_10_plot,Tleaf_12_plot,Tleaf_14_plot,Tleaf_16_plot, nrow=2)
ggsave("Tleaf_plot.png", plot = Tleaf_plot, dpi = 300, height = 5.0)


##


ggplot(dat_14_2,(aes(x = Construct2, y = PhiPS2, fill = Construct2)))+
  geom_boxplot()

###################
#summarize means by construct and hour
rd2_construct_df <- as.data.frame(aggregate(gsw ~ construct2 + time.hour., dat_rd2, mean))
write.csv(rd2_construct_df, file = "small_b_droughtstudy/round2/rd2_construct_means.csv")
#hour as factor

dat_rd2$time.hour.<- as.factor(dat_rd2$time.hour.)

rd_2_construct_plot <- ggplot(dat_rd2, aes(x = time.hour., y = gsw, fill = construct2))+
  geom_boxplot()+
  ylab("stomatal conductance (mol/m^2*sec)")
rd_2_construct_plot
ggsave(filename = "small_b_droughtstudy/round2/rd_2_construct_plot.png", plot = rd_2_construct_plot, dpi = 300, width = 10, height = 4, units = "in")


##predawn
library(janitor)

predawn_dat_rd2 <- read.csv("small_b_droughtstudy/round2/9_1_predawn_meas_proc.csv")
predawn_dat_rd2_li600 <- read.csv("small_b_droughtstudy/round2/predawn_sampling_9_1_22_proc.csv")
predawn_dat_rd2_li600 <- subset(predawn_dat_rd2_li600, select = c(ID,PhiPS2))

predawn_dat_rd2 <- inner_join(predawn_dat_rd2, predawn_dat_rd2_li600)
predawn_dat_rd2 <- inner_join(predawn_dat_rd2, dat_meta, by = "ID")

str(predawn_dat_rd2)
predawn_dat_rd2$MPa <- as.numeric(predawn_dat_rd2$MPa)
predawn_dat_rd2 <- na.omit(predawn_dat_rd2)

#heatmap with smoothing
#install.packages("latticeExtra")
library(latticeExtra)
str(predawn_dat_rd2)



#soil_moisture
rd_2_soil_moisture_plot <- levelplot(soil.moisture ~ row*column, predawn_dat_rd2, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_2_soil_moisture_plot

ggsave(filename = "small_b_droughtstudy/round2/rd_2_soilm_plot.png", plot = rd_2_soil_moisture_plot, dpi = 300)

#Water potential (MPa - some concerns over dew effects)
rd_2_mPA_plot <- levelplot(MPa ~ row*column, predawn_dat_rd2, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_2_mPA_plot
ggsave(filename = "small_b_droughtstudy/round2/rd_2_mPA_plot.png", plot = rd_2_mPA_plot, dpi = 300)

#PhiPS2
rd_2_fvfm_plot <- levelplot(PhiPS2 ~ row*column, predawn_dat_rd2, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))

rd_2_fvfm_plot
ggsave(filename = "small_b_droughtstudy/round2/rd_2_fvfm_plot.png", plot = rd_2_fvfm_plot, dpi = 300)

## correlation matrix
predawn_dat_rd2_corr <- subset(predawn_dat_rd2, select = c(soil.moisture,PhiPS2,MPa))
#predawn_dat_rd2_corr <- predawn_dat_rd2_corr %>% rename(Psileaf_MPa = MPa)
rd2_cor <-(cor(predawn_dat_rd2_corr, use = "complete.obs", method = "pearson"))

corrplot(rd2_cor, method = "number")



ggplot(predawn_dat_rd2,(aes(x = event_short, y = PhiPS2, fill = construct2)))+
  geom_boxplot()
