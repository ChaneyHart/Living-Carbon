####Analysis of diurnal curves on 6/27######



library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(tidyr)
library(stringi)

dat_6_27_6800 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/6_27_diurnal_datasheet_li6800_filled.csv")
dat_6_27_6800  <- subset(dat_6_27_6800, averaging == 15)

dat_6_15 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/dat_6_16_diurnal.csv")

dat_6_15_6800 <- subset(dat_6_15, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf))

dat_6_27_6800 <- subset(dat_6_27_6800, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf))
dat_6_27_6800_II <- dat_6_27_6800[c(1,2),]

#fix time stamp
stri_sub(dat_6_27_6800$date,1,8)
stri_sub(dat_6_27_6800$date,5,4) <- "/"
stri_sub(dat_6_27_6800$date,8,7) <- "/"
stri_sub(dat_6_27_6800_II$date,6,5) <- "20"
stri_sub(dat_6_15_6800$date,6,5) <- "20"
dat_6_27_6800$date
dat_6_27_6800_II$date
dat_6_15_6800$date
#substitute / for -
dat_6_27_6800$date <- gsub('\\/',"-",dat_6_27_6800$date)
dat_6_27_6800_II$date <- gsub('\\/',"-",dat_6_27_6800_II$date)
dat_6_15_6800$date <- gsub('\\/',"-",dat_6_15_6800$date)

#convert to posix
dat_6_27_6800$date <- ymd_hms(dat_6_27_6800$date)
dat_6_27_6800_II$date <- mdy_hm(dat_6_27_6800_II$date)
dat_6_15_6800$date <- mdy_hm(dat_6_15_6800$date)

#add info about time of sampling
dat_6_27_6800$hour <- hour(dat_6_27_6800$date)
dat_6_27_6800_II$hour <- hour(dat_6_27_6800_II$date)
dat_6_15_6800$hour <- hour(dat_6_15_6800$date)

#join datasets together
dat_6_27_6800 <- rbind(dat_6_27_6800,dat_6_27_6800_II)


#drop first two rows that lack proper dates
dat_6_27_6800 <- dat_6_27_6800[c(-1,-2),]

dat_6_27_6800 <- dat_6_27_6800 %>% mutate(Timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 15 & hour < 18 ~ 16,
  hour >= 18 & hour < 20 ~ 18,
  hour >= 20 & hour < 24 ~ 20
))

str(dat_6_27_6800)
dat_6_27_6800$PhiPS2 <- as.numeric(dat_6_27_6800$PhiPS2)
dat_6_27_6800$ETR <- as.numeric(dat_6_27_6800$ETR)
dat_6_27_6800$NPQ <- as.numeric(dat_6_27_6800$NPQ)
dat_6_27_6800$Fo. <- as.numeric(dat_6_27_6800$Fo.)
dat_6_27_6800$qP <- as.numeric(dat_6_27_6800$qP)
dat_6_27_6800$qN <- as.numeric(dat_6_27_6800$qN)
dat_6_27_6800$qL <- as.numeric(dat_6_27_6800$qL)
dat_6_27_6800$Timepoint <- as.factor(dat_6_27_6800$Timepoint)

dat_6_27_600 <- read.csv((file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/6_27_diurnal_li600 (filled).csv"))

dat_6_27_600 <- subset(dat_6_27_600, select = c(Time,Date,ID,Timepoint,gsw,VPDleaf,Fs,Fm.,PhiPS2,ETR,Qamb,Tleaf))
str(dat_6_27_600)
dat_6_27_600$Timepoint <- as.factor(dat_6_27_600$Timepoint)

#dat_6_27 <- full_join(dat_6_27_600, dat_6_6800, by = c("ID","Timepoint"))

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))

dat_6_27_6800 <- inner_join(dat_6_27_6800, LC_meta, by = "ID")
dat_6_27_600 <- inner_join(dat_6_27_600, LC_meta, by = "ID")
str(dat_6_27_600)
#add class factor

dat_6_27_6800 <- dat_6_27_6800 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

dat_6_27_600 <- dat_6_27_600 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control"))

#Create Summary table and graph



June_6800_summary_event <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(event_short,Timepoint) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
)

write.csv(June_6800_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_event.csv)")

June_6800_summary_tree <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(ID,Timepoint) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
)

write.csv(June_6800_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_tree.csv)")


June_600_summary_event <- subset(dat_6_27_600,) %>% dplyr::group_by(event_short,Timepoint) %>% dplyr::summarize(
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

write.csv(June_600_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_event.csv)")

June_600_summary_tree <- subset(dat_6_27_600,) %>% dplyr::group_by(ID,Timepoint) %>% dplyr::summarize(
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
  PhiPS2 = mean(PhiPS2)
)

write.csv(June_600_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_tree.csv)")


June_6800_summary_class <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(Class,Timepoint) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
)

write.csv(June_6800_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_class.csv)")


June_600_summary_class <- subset(dat_6_27_600,) %>% dplyr::group_by(Class,Timepoint) %>% dplyr::summarize(
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

write.csv(June_600_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_class.csv)")

# join datasets together

June_summary <- inner_join(June_600_summary, June_6800_summary, by = c("event_short","Timepoint"))

June_summary_class <- inner_join(June_600_summary_class, June_6800_summary_class, by = c("Class","Timepoint"))


#graph
#define colors
library(RColorBrewer)
my_colors <- c("gray0","indianred3","chartreuse4")

library(ggsignif)
library(gridExtra)
June_assimilation <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)

June_assimilation

June_gsw <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")
June_gsw

June_ETR <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")
June_ETR

June_diurnal_phys_plot <- grid.arrange(June_assimilation, June_gsw, June_ETR, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_diurnal_phys_plot.png", plot = June_diurnal_phys_plot, height = 10, units = "in", dpi = 300)

June_Tleaf <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = Tleaf), size = 2)+
  geom_line(aes(y = Tleaf))+ 
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")
June_Tleaf

June_VPDleaf <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = VPDleaf), size = 2)+
  geom_line(aes(y = VPDleaf))+ 
  geom_errorbar(aes(ymin = VPDleaf-VPDleaf_se, ymax = VPDleaf+VPDleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf VPD (kPa)")
June_VPDleaf

June_PAR <- ggplot(June_summary_class, aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PAR), size = 2)+
  geom_line(aes(y = PAR))+ 
  geom_errorbar(aes(ymin = PAR-PAR_se, ymax = PAR+PAR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*s)")
June_PAR

June_diurnal_conditions_plot <- grid.arrange(June_PAR, June_Tleaf, June_VPDleaf, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_diurnal_conditions_plot.png", plot = June_diurnal_conditions_plot, height = 10, units = "in", dpi = 300)

#read in soil moisture data

June_SM <- read.csv("LC_2023/2023_physiology_analysis/diurnal_experiments/6_27_moisture_data (filled).csv")
mean(June_SM$SM.Av)
mean(June_SM$SM.Av) - 2*(sd(June_SM$SM.Av)/(sqrt(48)))
mean(June_SM$SM.Av) + 2*(sd(June_SM$SM.Av)/(sqrt(48)))



