###Analysis of diurnal curves on 7/26######



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

dat_7_26_6800 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/Li6800_diurnal_7_26_23_filled.csv")
dat_7_26_6800  <- subset(dat_7_26_6800, averaging == 15)


dat_7_26_6800 <- subset(dat_7_26_6800, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf,A_dark))
#fix time stamp
stri_sub(dat_7_26_6800$date,1,8)
stri_sub(dat_7_26_6800$date,5,4) <- "/"
stri_sub(dat_7_26_6800$date,8,7) <- "/"
dat_7_26_6800$date
dat_7_26_6800$date <- gsub('\\/',"-",dat_7_26_6800$date)
dat_7_26_6800$date <- ymd_hms(dat_7_26_6800$date)

#add info about time of sampling
dat_7_26_6800$hour <- hour(dat_7_26_6800$date)


dat_7_26_6800 <- dat_7_26_6800 %>% mutate(Timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 16 & hour < 18 ~ 16,
  hour >= 18 & hour < 20 ~ 18,
  hour >= 20 & hour < 24 ~ 20
))

str(dat_7_26_6800)
dat_7_26_6800$PhiPS2 <- as.numeric(dat_7_26_6800$PhiPS2)
dat_7_26_6800$ETR <- as.numeric(dat_7_26_6800$ETR)
dat_7_26_6800$NPQ <- as.numeric(dat_7_26_6800$NPQ)
dat_7_26_6800$Fo. <- as.numeric(dat_7_26_6800$Fo.)
dat_7_26_6800$qP <- as.numeric(dat_7_26_6800$qP)
dat_7_26_6800$qN <- as.numeric(dat_7_26_6800$qN)
dat_7_26_6800$qL <- as.numeric(dat_7_26_6800$qL)
dat_7_26_6800$Timepoint <- as.factor(dat_7_26_6800$Timepoint)

dat_7_26_600 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/Li600_diurnal_7_26_23_filled.csv", header = FALSE)
dat_7_26_600 <- dat_7_26_600[c(-1,-3),]
dat_7_26_600 <- row_to_names(dat_7_26_600, 1)

dat_7_26_600 <- subset(dat_7_26_600, select = c(Time,Date,ID,gsw,VPDleaf,Fs,Fm,PhiPS2,ETR,Qamb,Tleaf))
str(dat_7_26_600)

#fix time format and add timepoint

dat_7_26_600$Time <- hms(dat_7_26_600$Time)
dat_7_26_600$hour <- hour(dat_7_26_600$Time)
dat_7_26_600 <- dat_7_26_600 %>% mutate(Timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 16 & hour < 18 ~ 16,
  hour >= 18 & hour < 20 ~ 18,
  hour >= 20 & hour < 24 ~ 20
))

dat_7_26_600$gsw <- as.numeric(dat_7_26_600$gsw)
dat_7_26_600$VPDleaf <- as.numeric(dat_7_26_600$VPDleaf)
dat_7_26_600$Fs <- as.numeric(dat_7_26_600$Fs)
dat_7_26_600$Fm <- as.numeric(dat_7_26_600$Fm)
dat_7_26_600$PhiPS2 <- as.numeric(dat_7_26_600$PhiPS2)
dat_7_26_600$ETR <- as.numeric(dat_7_26_600$ETR)
dat_7_26_600$Qamb <- as.numeric(dat_7_26_600$Qamb)
dat_7_26_600$Tleaf <- as.numeric(dat_7_26_600$Tleaf)
dat_7_26_600$Timepoint <- as.factor(dat_7_26_600$Timepoint)

#dat_7_26 <- left_join(dat_7_26_600, dat_7_26_6800, by = c("ID","Timepoint"))

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))

dat_7_26_6800 <- inner_join(dat_7_26_6800, LC_meta, by = "ID")
dat_7_26_600 <- inner_join(dat_7_26_600, LC_meta, by = "ID")

dat_7_26_6800 <- dat_7_26_6800 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

dat_7_26_600 <- dat_7_26_600 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

#define control trees and drought trees 
control_trees <- read.csv(file="LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
drought_trees <- LC_meta[!(LC_meta$ID %in% control_trees$ID),]


dat_7_26_6800 <- dat_7_26_6800 %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))


dat_7_26_600 <- dat_7_26_600 %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))


#How well do same measurements from different machines match
ggplot(dat_7_26, aes(x=gsw.x, y=gsw.y))+
  geom_point()+
  xlab("gsw_li600")+
  ylab("gsw_li6800")
ggplot(dat_7_26, aes(x=VPDleaf.x, y=VPDleaf.y))+
  geom_point()             
ggplot(dat_7_26, aes(x=PhiPS2.x, y=PhiPS2.y))+
  geom_point()  
ggplot(dat_7_26, aes(x=ETR.x, y=ETR.y))+
  geom_point()  


#Explore event differences and hourly trends

ggplot(dat_7_26,aes(x= Timepoint, y=A, fill = construct2))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=A, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=gsw.x, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=gsw.y, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=PhiPS2.x, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=PhiPS2.y, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=ETR.x, fill = Class))+
  geom_boxplot()

ggplot(dat_7_26,aes(x= Timepoint, y=ETR.y, fill = Class))+
  geom_boxplot()

###
# Summary tables and graphs

July_6800_summary_event <- subset(dat_7_26_6800, A > 0) %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
)

write.csv(July_6800_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_event.csv)")

July_6800_summary_tree <- subset(dat_7_26_6800, A > 0) %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  Fs = mean(Fs),
  Fm_prime = mean(Fm.),
  Fo. = mean(Fo.),
  A_dark = mean(A_dark),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  PAR_sd = sd(Qin),
  PAR = mean(Qin),
  PAR_se = (PAR_sd/(sqrt(n))),
  PhiPS2 = mean(PhiPS2)
  
  
)

write.csv(July_6800_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_tree.csv)")


July_600_summary_event <- subset(dat_7_26_600,) %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n))),
  PAR_sd = sd(Qamb),
  PAR = mean(Qamb),
  PAR_se = (PAR_sd/(sqrt(n))),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  VPDleaf_sd = sd(VPDleaf),
  VPDleaf = mean(VPDleaf),
  VPDleaf_se = (VPDleaf_sd/(sqrt(n))),
)

write.csv(July_600_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_event.csv)")

July_600_summary_tree <- subset(dat_7_26_600,) %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n))),
  PAR_600_sd = sd(Qamb),
  PAR_600 = mean(Qamb),
  PAR_600_se = (PAR_600_sd/(sqrt(n))),
  ETR_600_sd = sd(ETR),
  ETR_600 = mean(ETR),
  ETR_600_se = (ETR_600_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  VPDleaf_sd = sd(VPDleaf),
  VPDleaf = mean(VPDleaf),
  VPDleaf_se = (VPDleaf_sd/(sqrt(n))),
  PhiPS2_600 = mean(PhiPS2)
)

write.csv(July_600_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_tree.csv)")

July_6800_summary_class <- subset(dat_7_26_6800, A > 0) %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
)

write.csv(July_6800_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_class.csv)")


July_600_summary_class <- subset(dat_7_26_600,) %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n))),
  PAR_sd = sd(Qamb),
  PAR = mean(Qamb),
  PAR_se = (PAR_sd/(sqrt(n))),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  VPDleaf_sd = sd(VPDleaf),
  VPDleaf = mean(VPDleaf),
  VPDleaf_se = (VPDleaf_sd/(sqrt(n))),
)

write.csv(July_600_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_class.csv)")


# join datasets together

July_summary <- inner_join(July_600_summary, July_6800_summary, by = c("event_short","Timepoint"))

July_summary_class <- inner_join(July_600_summary_class, July_6800_summary_class, by = c("Class","Timepoint","treatment"))


#graph
#define colors
library(RColorBrewer)
my_colors <- c("gray0","indianred3","chartreuse4")

library(ggsignif)
library(gridExtra)
July_assimilation <- ggplot(July_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)

July_assimilation

July_gsw <- ggplot(July_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")
July_gsw

July_ETR <- ggplot(July_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")
July_ETR

July_diurnal_phys_plot <- grid.arrange(July_assimilation, July_gsw, July_ETR, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_diurnal_phys_plot.png", plot = July_diurnal_phys_plot, height = 10, units = "in", dpi = 300)

July_Tleaf <- ggplot(July_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = Tleaf), size = 2)+
  geom_line(aes(y = Tleaf))+ 
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")
July_Tleaf

July_VPDleaf <- ggplot(July_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = VPDleaf), size = 2)+
  geom_line(aes(y = VPDleaf))+ 
  geom_errorbar(aes(ymin = VPDleaf-VPDleaf_se, ymax = VPDleaf+VPDleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf VPD (kPa)")
July_VPDleaf

July_PAR <- ggplot(July_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PAR), size = 2)+
  geom_line(aes(y = PAR))+ 
  geom_errorbar(aes(ymin = PAR-PAR_se, ymax = PAR+PAR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*s)")
July_PAR

July_diurnal_conditions_plot <- grid.arrange(July_PAR, July_Tleaf, July_VPDleaf, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_diurnal_conditions_plot.png", plot = July_diurnal_conditions_plot, height = 10, units = "in", dpi = 300)



##Examine subset with full fluorescent data

July_fluor_set <- inner_join(July_6800_summary_tree,July_600_summary_tree, by = c("ID","treatment","Timepoint"))

#define respiration_data

July_resp_dat <- subset(dat_7_26_6800, A < 0)
July_resp_dat$Rd <- July_resp_dat$A_dark
July_resp_dat$Resp_temp <- July_resp_dat$Tleaf

July_resp_dat <- subset(July_resp_dat, select = c(ID,Rd,Resp_temp))
July_fluor_set <- left_join(July_fluor_set, July_resp_dat)

July_fluor_set$R <- (July_fluor_set$Rd*2.2)^((July_fluor_set$Tleaf - July_fluor_set$Resp_temp)/10)

#(July_fluor_set$Rd*2.2)^July_fluor_set$resp_factor
#July_fluor_set$resp_factor <- (July_fluor_set$Tleaf - July_fluor_set$Resp_temp)/10

July_fluor_set$Theta_e <- (4*(July_fluor_set$PhiPS2_600 - 0.026947))/7.567127
July_fluor_set$Jt <- July_fluor_set$Theta_e*July_fluor_set$PAR_600
July_fluor_set$Jo <- (2/3)*(July_fluor_set$Jt - (4*(July_fluor_set$A + July_fluor_set$Rd)))
July_fluor_set$Rl <- (1/12)*(July_fluor_set$Jt - (4*(July_fluor_set$A + July_fluor_set$Rd)))
July_fluor_set$Jc <- (1/3)*(July_fluor_set$Jt + (8*(July_fluor_set$A + July_fluor_set$Rd)))
July_fluor_set$Jo_ratio <- July_fluor_set$Jo/July_fluor_set$Jt

July_fluor_set <- inner_join(July_fluor_set, LC_meta, by = "ID")

July_fluor_set <- July_fluor_set %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

July_fluor_set <- July_fluor_set %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))

July_fluor_set <- subset(July_fluor_set, PAR_600 > 300)

July_fluor_set_summary_event <- July_fluor_set %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  PAR_sd = sd(PAR),
  PAR = mean(PAR),
  PAR_se = (PAR_sd/(sqrt(n))),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n))),
  Jt_sd = sd(Jt),
  Jt = mean(Jt),
  Jt_se = (Jt_sd/(sqrt(n))),
  Jo_sd = sd(Jo),
  Jo = mean(Jo),
  Jo_se = (Jo_sd/(sqrt(n))),
  Jc_sd = sd(Jc),
  Jc = mean(Jc),
  Jc_se = (Jc_sd/(sqrt(n))),
  Rd_sd = sd(Rd),
  Rd = mean(Rd),
  Rd_se = (Rd_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  Jo_ratio_sd = sd(Jo_ratio),
  Jo_ratio = mean(Jo_ratio),
  Jo_ratio_se = (Jo_ratio_sd/(sqrt(n))),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n)))
)  



write.csv(July_fluor_set_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_event.csv)")

July_fluor_set_summary_class <- July_fluor_set %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  PAR_sd = sd(PAR),
  PAR = mean(PAR),
  PAR_se = (PAR_sd/(sqrt(n))),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n))),
  Jt_sd = sd(Jt),
  Jt = mean(Jt),
  Jt_se = (Jt_sd/(sqrt(n))),
  Jo_sd = sd(Jo),
  Jo = mean(Jo),
  Jo_se = (Jo_sd/(sqrt(n))),
  Jc_sd = sd(Jc),
  Jc = mean(Jc),
  Jc_se = (Jc_sd/(sqrt(n))),
  Rd_sd = sd(Rd),
  Rd = mean(Rd),
  Rd_se = (Rd_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  Jo_ratio_sd = sd(Jo_ratio),
  Jo_ratio = mean(Jo_ratio),
  Jo_ratio_se = (Jo_ratio_sd/(sqrt(n))),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n)))
)

write.csv(July_fluor_set_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_class.csv)")

July_fluor_set_summary_tree <- July_fluor_set %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = (ETR_sd/(sqrt(n))),
  PAR_sd = sd(PAR),
  PAR = mean(PAR),
  PAR_se = (PAR_sd/(sqrt(n))),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n))),
  Jt_sd = sd(Jt),
  Jt = mean(Jt),
  Jt_se = (Jt_sd/(sqrt(n))),
  Jo_sd = sd(Jo),
  Jo = mean(Jo),
  Jo_se = (Jo_sd/(sqrt(n))),
  Jc_sd = sd(Jc),
  Jc = mean(Jc),
  Jc_se = (Jc_sd/(sqrt(n))),
  Rd_sd = sd(Rd),
  Rd = mean(Rd),
  Rd_se = (Rd_sd/(sqrt(n))),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = (Tleaf_sd/(sqrt(n))),
  Jo_ratio_sd = sd(Jo_ratio),
  Jo_ratio = mean(Jo_ratio),
  Jo_ratio_se = (Jo_ratio_sd/(sqrt(n))),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = (gsw_sd/(sqrt(n)))
)

write.csv(July_fluor_set_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_tree.csv)")

