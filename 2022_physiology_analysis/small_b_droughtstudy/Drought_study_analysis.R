library(ggplot2)
library(lubridate)
library(gridExtra)
library((ggthemes))
library(ggsignif)
library(janitor)
library(dplyr)
library(corrplot)

#import data
#####################################
#round 1

#Diurnal curves (8_15_22)
dat_rd1_Still <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round1/8_15_drought_00299_proc_chaneys_li600.csv", header = TRUE)
dat_rd1_Still <- subset(dat_rd1_Still, select = c(Time,Hour,ID,gsw,VPDleaf,PhiPS2,ETR,Tleaf,Qamb))
dat_rd1_LC <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round1/8_15_drought_00354_proc_janes_li600.csv", header = TRUE)
dat_rd1_LC <- subset(dat_rd1_LC, select = c(Time,Hour,ID,gsw,VPDleaf,PhiPS2,ETR,Tleaf,Qamb))
#combine data from both instruments
dat_rd1_combined <- rbind(dat_rd1_Still, dat_rd1_LC)
#import meta data w/ event, construct...etc
dat_meta <- read.csv("2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/DBH_H_timeline_CT1_excluded_9_22.csv")
dat_meta <- subset(dat_meta, select = c(row,column,ID,event_short, construct, construct2,block,H419))
#match id to metadata
dat_rd1 <- inner_join(dat_rd1_combined, dat_meta, by = "ID")
dat_rd1$Hour <- as.factor(dat_rd1$Hour)


#summarize means by construct and hour
rd1_construct_df <- as.data.frame(aggregate(gsw ~ construct2 + Hour, dat_rd1, mean))
write.csv(rd1_construct_df, file = "small_b_droughtstudy/round1/rd1_construct_means.csv")
#graph
rd_1_construct_plot <- ggplot(dat_rd1, aes(x = Hour, y = gsw, fill = construct2))+
  geom_boxplot()+
  ylab("stomatal conductance (mol/m^2*sec)")

rd_1_construct_plot
ggsave(filename = "small_b_droughtstudy/round1/rd_1_construct_plot.png", plot = rd_1_construct_plot, dpi = 300, width = 10, height = 4, units = "in")

#summarize means by event and hour
rd1_event_df <- as.data.frame(aggregate(gsw ~ event_short + Hour, dat_rd1, mean))
ggplot(dat_rd1, aes(x = Hour, y = gsw, fill = event_short))+
  geom_boxplot()

###Look at each hour seperately 
#10###3
dat_10 <- subset(dat_rd1, Hour =="10")

gsw_10_plot <- ggplot(dat_10, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (µmol/mol s)")+
  ggtitle("10:00")+
  theme(legend.position = "None")

gsw_10_plot2 <- ggplot(dat_10, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("10:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.8)+
  theme(legend.position = "None")

gsw_10_plot <- gsw_10_plot + colorScale
gsw_10_plot
gsw_10_plot2 
ggsave("small_b_droughtstudy/round1/gsw_10_rd1.png", plot = gsw_10_plot2,dpi = 300)

####
lm_10_rd1 <- lm(gsw~construct2, dat_10)
anova(lm_10_rd1)



####12#####
dat_12 <- subset(dat_rd1, Hour =="12")

gsw_12_plot <- ggplot(dat_12, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("12:00")+
  theme(legend.position = "None")

gsw_12_plot2 <- ggplot(dat_12, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("12:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.8)+
  theme(legend.position = "None")

gsw_12_plot <- gsw_12_plot + colorScale
gsw_12_plot
gsw_12_plot2 
ggsave("small_b_droughtstudy/round1/gsw_12_rd1.png", plot = gsw_12_plot2,dpi = 300)

###
lm_12_rd1 <- lm(gsw~construct2, dat_12)
anova(lm_12_rd1)

######14############

dat_14 <- subset(dat_rd1, Hour =="14")

gsw_14_plot <- ggplot(dat_14, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (µmol/mol s)")+
  ggtitle("14:00")+
  theme(legend.position = "None")

gsw_14_plot2 <- ggplot(dat_14, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("14:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.65)+
  theme(legend.position = "None")

gsw_14_plot <- gsw_14_plot + colorScale
gsw_14_plot
gsw_14_plot2 
ggsave("small_b_droughtstudy/round1/gsw_14_rd1.png", plot = gsw_14_plot2,dpi = 300)

#stats
lm_14_rd1 <- lm(gsw~construct2, dat_14)
anova(lm_14_rd1)

######16###########

dat_16 <- subset(dat_rd1, Hour =="16")

gsw_16_plot <- ggplot(dat_16, aes(x=reorder(event_short, gsw),y = gsw, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("gsw (µmol/mol s)")+
  ggtitle("16:00")+
  theme(legend.position = "None")

gsw_16_plot2 <- ggplot(dat_16, aes(x=reorder(construct2, gsw),y = gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  ggtitle("16:00")+
  theme(legend.position = "None")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.8)+
  theme(legend.position = "None")

gsw_16_plot <- gsw_16_plot + colorScale
gsw_16_plot
gsw_16_plot2 
ggsave("small_b_droughtstudy/round1/gsw_16_rd1.png", plot = gsw_16_plot2,dpi = 300)

rd1_gsw_plot2 <- grid.arrange(gsw_10_plot2,gsw_12_plot2, gsw_14_plot2,gsw_16_plot2, nrow=2)
ggsave("small_b_droughtstudy/round1/gsw_plot2_rd1.png", plot = rd1_gsw_plot2,height = 5.0, dpi = 300)

####
lm_16_rd1 <- lm(gsw~construct2, dat_16)
anova(lm_16_rd1)


#############
#stats
rd1_gsw_anova <- lm(gsw ~ construct2 + Hour, dat_rd1)
summary(rd1_gsw_anova)
anova(rd1_gsw_anova)


########predawn sampling##########
#import data
predawn_dat_rd1 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round1/8_19_predawn_meas.csv")
predawn_dat_rd1_li600 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round1/Night sampling/Night_sampling_CH_2022_08_19_proc.csv")
predawn_dat_rd1_li600 <- subset(predawn_dat_rd1_li600, select = c(ID,Fv.Fm))

predawn_dat_rd1 <- inner_join(predawn_dat_rd1, predawn_dat_rd1_li600)
predawn_dat_rd1 <- inner_join(predawn_dat_rd1, dat_meta, by = "ID")

predawn_dat_rd1$soil.moisture <- as.numeric(predawn_dat_rd1$soil.moisture)
predawn_dat_rd1$Fv.Fm <- as.numeric(predawn_dat_rd1$Fv.Fm)
predawn_dat_rd1$MPa <- as.numeric(predawn_dat_rd1$MPa)

predawn_dat_rd1 <- na.omit(predawn_dat_rd1)
#heat map of soil moisture
#rough first attempt
ggplot(predawn_dat_rd1, aes(x = column, y =row))+
  geom_point(aes(color=soil.moisture))

ggplot(predawn_dat_rd1, aes(row, column, fill = soil.moisture))+
  geom_density_2d_filled()+
  stat_density2d(geom = "tile", aes(fill =..density..), contour = FALSE)

#with smoothing
#install.packages("latticeExtra")
library(latticeExtra)
str(predawn_dat_rd1)


#soil moisture

rd_1_soil_moisture_plot <- levelplot(soil.moisture ~ row*column, predawn_dat_rd1, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_1_soil_moisture_plot
ggsave(filename = "small_b_droughtstudy/round1/rd_1_soilm_plot.png", plot = rd_1_soil_moisture_plot, dpi = 300)

#Predawn water potential
rd_1_mPA_plot <- levelplot(MPa ~ row*column, predawn_dat_rd1, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_1_mPA_plot
ggsave(filename = "small_b_droughtstudy/round1/rd_1_mPA_plot.png", plot = rd_1_mPA_plot, dpi = 300)

#Fv/Fm
rd_1_fvfm_plot <- levelplot(Fv.Fm ~ row*column, predawn_dat_rd1, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_1_fvfm_plot
ggsave(filename = "small_b_droughtstudy/round1/rd_1_fvfm_plot.png", plot = rd_1_fvfm_plot, dpi = 300)

## correlation matrix
predawn_dat_rd1_corr <- subset(predawn_dat_rd1, select = c(soil.moisture,Fv.Fm,MPa))

rd1_cor <-(cor(predawn_dat_rd1_corr, use = "complete.obs", method = "pearson"))
corrplot(rd1_cor, method = "number")


ggplot(predawn_dat_rd1,(aes(x = reorder(event_short,Fv.Fm), y = Fv.Fm, fill = construct2)))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Fv/Fm")

rd_1_fvfm_model <- lm(Fv.Fm ~ event_short, predawn_dat_rd1)
anova(rd_1_fvfm_model)


#####
#create plots that map each variable over the whole month

# read in data for rounds 2 & 3

predawn_dat_rd3 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round3/rd_3_predawn_meas_proc.csv")
predawn_dat_rd3_li600 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round3/9_22_predawn_li600_rd3_proc.csv")
predawn_dat_rd3_li600 <- subset(predawn_dat_rd3_li600, select = c(ID,Fv.Fm))
predawn_dat_soil_moisture <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round3/9_19_soilmoisture.csv")

predawn_dat_rd3 <- inner_join(predawn_dat_rd3, predawn_dat_rd3_li600)
predawn_dat_rd3 <- left_join(predawn_dat_rd3, predawn_dat_soil_moisture, by = "ID")
predawn_dat_rd3 <- left_join(predawn_dat_rd3, dat_meta, by = "ID")
predawn_dat_rd3 <- subset(predawn_dat_rd3, select = c(ID,obs,MPa,bar,Fv.Fm,soil.moisture,row,column,event_short,construct,construct2,H419))

predawn_dat_rd3 <- na.omit(predawn_dat_rd3)

###

predawn_dat_rd2 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round2/9_1_predawn_meas_proc.csv")
predawn_dat_rd2_li600 <- read.csv("2022_physiology_analysis/small_b_droughtstudy/round2/predawn_sampling_9_1_22_proc.csv")
predawn_dat_rd2_li600 <- subset(predawn_dat_rd2_li600, select = c(ID,PhiPS2))

predawn_dat_rd2 <- inner_join(predawn_dat_rd2, predawn_dat_rd2_li600)
predawn_dat_rd2 <- inner_join(predawn_dat_rd2, dat_meta, by = "ID")

str(predawn_dat_rd2)
predawn_dat_rd2$MPa <- as.numeric(predawn_dat_rd2$MPa)
predawn_dat_rd2 <- na.omit(predawn_dat_rd2)


###
#summarize event means

predawn_rd1_construct <- predawn_dat_rd1 %>% group_by(construct2)

predawn_rd1_summary <- predawn_rd1_construct %>% dplyr::summarise(
  n = n(),
  rd1_SM_sd = sd(soil.moisture),
  rd1_SM = mean(soil.moisture),
  rd1_SM_se = (rd1_SM_sd/(sqrt(n))),
  rd1_MPa_sd = sd(MPa),
  rd1_MPa = mean(MPa),
  rd1_MPa_se = (rd1_MPa_sd/(sqrt(n))),
  rd1_FvFm_sd = sd(Fv.Fm),
  rd1_FvFm = mean(Fv.Fm),
  rd1_FvFm_se = (rd1_FvFm_sd/(sqrt(n)))
)

predawn_rd2_construct <- predawn_dat_rd2 %>% group_by(construct2)

predawn_rd2_summary <- predawn_rd2_construct %>% dplyr::summarise(
  n = n(),
  rd2_SM_sd = sd(soil.moisture),
  rd2_SM = mean(soil.moisture),
  rd2_SM_se = (rd2_SM_sd/(sqrt(n))),
  rd2_MPa_sd = sd(MPa),
  rd2_MPa = mean(MPa),
  rd2_MPa_se = (rd2_MPa_sd/(sqrt(n))),
  rd2_FvFm_sd = sd(fv.fm),
  rd2_FvFm = mean(fv.fm),
  rd2_FvFm_se = (rd2_FvFm_sd/(sqrt(n)))
)

predawn_rd3_construct <- predawn_dat_rd3 %>% group_by(construct2)

predawn_rd3_summary <- predawn_rd3_construct %>% dplyr::summarise(
  n = n(),
  rd3_SM_sd = sd(soil.moisture),
  rd3_SM = mean(soil.moisture),
  rd3_SM_se = (rd3_SM_sd/(sqrt(n))),
  rd3_MPa_sd = sd(MPa),
  rd3_MPa = mean(MPa),
  rd3_MPa_se = (rd3_MPa_sd/(sqrt(n))),
  rd3_FvFm_sd = sd(Fv.Fm),
  rd3_FvFm = mean(Fv.Fm),
  rd3_FvFm_se = (rd3_FvFm_sd/(sqrt(n)))
)

#xport to compile in excel
write.csv(predawn_rd1_summary, file = "2022_physiology_analysis/small_b_droughtstudy/predawn_rd1_summary.csv")
write.csv(predawn_rd2_summary, file = "2022_physiology_analysis/small_b_droughtstudy/predawn_rd2_summary.csv")
write.csv(predawn_rd3_summary, file = "2022_physiology_analysis/small_b_droughtstudy/predawn_rd3_summary.csv")

#read compiled file back in
predawn_summary <- read.csv("2022_physiology_analysis/small_b_droughtstudy/Drought_study_conditions_summary.csv")
predawn_summary$Days <- as.factor(predawn_summary$Days)

SM_plot <- ggplot(predawn_summary, aes(x=Days,y=Soil.Moisture))+
  geom_point(aes(color = construct),size = 4)+
  geom_errorbar(aes(ymin=Soil.Moisture-SM_se, ymax=Soil.Moisture+SM_se),linewidth = 0.4, width = 0.3)+
  ylab("Soil moisture (m^3/m^3)")+
  xlab("Days since drought")

SM_plot
ggsave(filename = "SM.png", plot=SM_plot, width=4.5, height =4, units = "in",dpi = 300)


MPa_plot <- ggplot(predawn_summary, aes(x=Days,y=Mpa))+
  geom_point(aes(color = construct),size = 4)+
  geom_errorbar(aes(ymin=Mpa-Mpa_se, ymax=Mpa+Mpa_se),linewidth = 0.4, width = 0.3)+
  ylab("Leaf water potential (MPa)")+
  xlab("Days since drought")
MPa_plot

ggsave(filename = "MPa.png", plot=MPa_plot, width=4.5, height = 4,units = "in",dpi = 300)

FvFm_plot <- ggplot(predawn_summary, aes(x=Days,y=FvFm))+
  geom_point(aes(color = construct),size = 4)+
  geom_errorbar(aes(ymin=FvFm-FvFm_se, ymax=FvFm+FvFm_se),linewidth = 0.4, width = 0.3)+
  ylab("Max quantum yield of PSII (Fv/Fm)")+
  xlab("Days since drought")
FvFm_plot

#Remove day 17 data

predawn_summary_II <- subset(predawn_summary, Days != 17)

FvFm_plot_II <- ggplot(predawn_summary_II, aes(x=Days,y=FvFm))+
  geom_point(aes(color = construct),size = 4)+
  geom_errorbar(aes(ymin=FvFm-FvFm_se, ymax=FvFm+FvFm_se),linewidth = 0.4, width = 0.3)+
  ylab("Max quantum yield of PSII (Fv/Fm)")+
  xlab("Days since drought")
FvFm_plot_II


ggsave(filename = "FvFm.png", plot=FvFm_plot_II, width=4.5, height = 4, units = "in",dpi = 300)
       