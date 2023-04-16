library(ggplot2)
library(lubridate)
library(gridExtra)
library((ggthemes))
library(ggsignif)
library(janitor)
library(dplyr)

#import data
#####################################
#round 1

#Diurnal curves (8_15_22)
dat_rd3 <- read.csv("small_b_droughtstudy/round3/SB_rd3_li600/9_19_sb_600_proc.csv", header = TRUE)
dat_rd3 <- subset(dat_rd3, select = c(Time,ID,gsw,VPDleaf,PhiPS2,ETR,Tleaf,Qamb))

#import meta data w/ event, construct...etc
dat_meta <- read.csv("DBH_H_timeline_CT1_excluded_9_22.csv")
dat_meta <- subset(dat_meta, select = c(row,column,ID,event_short, construct, construct2,block,H419))
#match id to metadata
dat_rd3 <- inner_join(dat_rd3, dat_meta, by = "ID")



#summarize means by construct 
rd3_construct_df <- as.data.frame(aggregate(gsw ~ construct2, dat_rd3, mean))
write.csv(rd3_construct_df, file = "small_b_droughtstudy/round3/rd3_construct_means.csv")
#graph
rd_3_construct_plot <- ggplot(dat_rd3, aes(x = construct2, y = gsw, fill = construct2))+
  geom_boxplot()+
  ylab("stomatal conductance (mol/m^2*sec)")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.5)+
  theme(legend.position = "None")

rd_3_construct_plot
ggsave(filename = "small_b_droughtstudy/round3/rd_3_construct_plot.png", plot = rd_3_construct_plot, width = 3, height = 3, units = "in",dpi = 300)

lm_rd3_gsw_12 <- lm(gsw~construct2, dat_rd3)
anova(lm_rd3_gsw_12)

#summarize means by event and hour
rd3_event_df <- as.data.frame(aggregate(gsw ~ event_short + Hour, dat_rd3, mean))
ggplot(dat_rd3, aes(x = Hour, y = gsw, fill = event_short))+
  geom_boxplot()

#stats
rd3_gsw_anova <- lm(gsw ~ construct2, dat_rd3)
summary(rd3_gsw_anova)
anova(rd3_gsw_anova)


########predawn sampling##########
#import data
predawn_dat_rd3 <- read.csv("small_b_droughtstudy/round3/rd_3_predawn_meas_proc.csv")
predawn_dat_rd3_li600 <- read.csv("small_b_droughtstudy/round3/9_22_predawn_li600_rd3_proc.csv")
predawn_dat_rd3_li600 <- subset(predawn_dat_rd3_li600, select = c(ID,Fv.Fm))
predawn_dat_soil_moisture <- read.csv("small_b_droughtstudy/round3/9_19_soilmoisture.csv")

predawn_dat_rd3 <- inner_join(predawn_dat_rd3, predawn_dat_rd3_li600)
predawn_dat_rd3 <- left_join(predawn_dat_rd3, predawn_dat_soil_moisture, by = "ID")
predawn_dat_rd3 <- left_join(predawn_dat_rd3, dat_meta, by = "ID")
predawn_dat_rd3 <- subset(predawn_dat_rd3, select = c(ID,obs,MPa,bar,Fv.Fm,soil.moisture,row,column,event_short,construct,construct2,H419))

predawn_dat_rd3 <- na.omit(predawn_dat_rd3)
#heat map of soil moisture
#rough first attempt
ggplot(predawn_dat_rd3, aes(row, column, fill = soil.moisture))+
  geom_tile()

ggplot(predawn_dat_rd3, aes(row, column, fill = soil.moisture))+
  geom_density_2d_filled()

#with smoothing
#install.packages("latticeExtra")
library(latticeExtra)
str(predawn_dat_rd3)

predawn_dat_rd3$soil.moisture <- as.numeric(predawn_dat_rd3$soil.moisture)
predawn_dat_rd3$Fv.Fm <- as.numeric(predawn_dat_rd3$Fv.Fm)
predawn_dat_rd3$mPA <- as.numeric(predawn_dat_rd3$mPA)



#soil moisture
rd_3_soil_moisture_plot <- levelplot(soil.moisture ~ row*column, predawn_dat_rd3, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_3_soil_moisture_plot
ggsave(filename = "small_b_droughtstudy/round3/rd_3_soilm_plot.png", plot = rd_1_soil_moisture_plot, dpi = 300)

#Predawn water potential
rd_3_mPA_plot <- levelplot(MPa ~ row*column, predawn_dat_rd3, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_3_mPA_plot
ggsave(filename = "small_b_droughtstudy/round3/rd_3_mPA_plot.png", plot = rd_1_mPA_plot, dpi = 300)

#Fv/Fm
rd_3_fvfm_plot <- levelplot(Fv.Fm ~ row*column, predawn_dat_rd3, panel = panel.levelplot.points, cex = 1.2)+
  layer_(panel.2dsmoother(..., n =200))
rd_3_fvfm_plot
ggsave(filename = "small_b_droughtstudy/round3/rd_3_fvfm_plot.png", plot = rd_1_fvfm_plot, dpi = 300)

## correlation matrix
predawn_dat_rd3_corr <- subset(predawn_dat_rd3, select = c(soil.moisture,Fv.Fm,MPa))

rd3_cor_df <- cor(predawn_dat_rd3_corr, use = "complete.obs", method = "pearson")
corrplot(rd3_cor_df, method = "number")


rd3_fvfm_construct_plot <- ggplot(predawn_dat_rd3,(aes(x = construct2, y = Fv.Fm, fill = construct2)))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.85)+
  theme(legend.position = "None")

rd3_fvfm_construct_plot

rd3_SM_construct_plot <- ggplot(predawn_dat_rd3,(aes(x = construct2, y = soil.moisture, fill = construct2)))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.2)+
  theme(legend.position = "None")

rd3_SM_construct_plot

rd3_mPA_construct_plot <- ggplot(predawn_dat_rd3,(aes(x = construct2, y = MPa, fill = construct2)))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = -0.3)+
  theme(legend.position = "None")

rd3_mPA_construct_plot

ggsave(filename = "small_b_droughtstudy/round3/rd_3_fvfm_construct_plot.png", plot = rd3_fvfm_construct_plot, dpi = 300)

#########assimilation data #####################

assimilation_dat <- read.csv("small_b_droughtstudy/round3/Drought_assimilation_rd3_proc.csv")
assimilation_dat <- subset(assimilation_dat, select = c(ID,A,Ci,gsw,VPDleaf,PhiPS2,ETR,PhiCO2,Tleaf,Qamb_out))
assimilation_dat <- inner_join(assimilation_dat, dat_meta, by="ID")

#assimilation
rd3_assim_plot <- ggplot(assimilation_dat, aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 23)+
  theme(legend.position = "None")

ggsave(filename = "small_b_droughtstudy/round3/rd3_assim_plot.png", plot = rd3_assim_plot, width = 3, height = 3, units = "in", dpi=300)


rd3_A_gsw <- ggplot(assimilation_dat, aes(x=gsw, y=A))+
  geom_point(aes(shape=construct2, color=construct2))+
  xlab("stomatal conductance (mol/m^2*s)")+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_smooth(method = "lm")

ggsave(filename = "rd3_A_gsw.png", plot = rd3_A_gsw, width = 4, height = 4, units = "in", dpi = 300)  

A_gsw_model <- lm(A~gsw, assimilation_dat)
summary(A_gsw_model)

rd3_gsw2_plot <- ggplot(assimilation_dat, aes(x=construct2, y=gsw, fill = construct2))+
  geom_boxplot()+
  xlab("Construct")+
  ylab("gsw (mol/m^2*s)")+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.4)+
  theme(legend.position = "None")

rd3_gsw2_plot
ggsave(filename = "rd3_gsw_9_22.png", plot = rd3_gsw2_plot, width = 3, height = 3, units = "in", dpi = 300)
# stats


ht_A <- lm(A ~ construct2, assimilation_dat)
anova(anova_drought_A)

anova_droughtA_event <- lm(A ~ event_short, assimilation_dat)
anova(anova_droughtA_event)

#phiPS2
ggplot(assimilation_dat, aes(x=construct2, y=PhiPS2, fill = construct2))+
  geom_boxplot()

ggplot(assimilation_dat, aes(x=construct2, y=gsw, fill = construct2))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("control","transgenic")),map_signif_level=TRUE, y_position = 0.4)+
  theme(legend.position = "None")

#ETR
ggplot(assimilation_dat, aes(x=event_short, y=ETR, fill = construct2))+
  geom_boxplot()

ggplot(assimilation_dat, aes(x=construct2, y=ETR, fill = construct2))+
  geom_boxplot()

#PhiCO2 to PhiPS2
ggplot(assimilation_dat, aes(x=PhiCO2, y =PhiPS2, color = event_short))+
  geom_point()

#combine Fv/Fm w/ phiPs2 to look at quenching
quenching_dat_1 <- read.csv("small_b_droughtstudy/round3/Drought_assimilation_rd3_proc.csv")
quenching_dat_1 <- subset(quenching_dat_1, select = c(ID,A,Ci,gsw,Fs,Fm.,PhiPS2,ETR,PhiCO2,Fo.,Qin,Tleaf))
quenching_dat_2 <- read.csv("small_b_droughtstudy/round3/9_22_predawn_li600_rd3_proc.csv")
quenching_dat_2 <- subset(quenching_dat_2, select = c(ID,Fv.Fm,Fo,Fm))

quenching_dat <- inner_join(quenching_dat_1, quenching_dat_2, by ="ID")
NPQ <- quenching_dat$Fm - quenching
