#importing and processing August diurnal data

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

dat_8_28_6800 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/diurnal_8_28_filled.csv")
dat_8_28_6800  <- subset(dat_8_28_6800, averaging == 5)

dat_8_28_6800 <- subset(dat_8_28_6800, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf,A_dark))
#fix time stamp
stri_sub(dat_8_28_6800$date,1,8)
stri_sub(dat_8_28_6800$date,5,4) <- "/"
stri_sub(dat_8_28_6800$date,8,7) <- "/"
dat_8_28_6800$date
dat_8_28_6800$date <- gsub('\\/',"-",dat_8_28_6800$date)
dat_8_28_6800$date <- ymd_hms(dat_8_28_6800$date)

#add info about time of sampling
dat_8_28_6800$hour <- hour(dat_8_28_6800$date)


dat_8_28_6800 <- dat_8_28_6800 %>% mutate(Timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 16 & hour < 18 ~ 16,
  hour >= 18 & hour < 20 ~ 18,
  hour >= 20 & hour < 24 ~ 20
))

str(dat_8_28_6800)
dat_8_28_6800$PhiPS2 <- as.numeric(dat_8_28_6800$PhiPS2)
dat_8_28_6800$ETR <- as.numeric(dat_8_28_6800$ETR)
dat_8_28_6800$NPQ <- as.numeric(dat_8_28_6800$NPQ)
dat_8_28_6800$Fo. <- as.numeric(dat_8_28_6800$Fo.)
dat_8_28_6800$qP <- as.numeric(dat_8_28_6800$qP)
dat_8_28_6800$qN <- as.numeric(dat_8_28_6800$qN)
dat_8_28_6800$qL <- as.numeric(dat_8_28_6800$qL)
dat_8_28_6800$Timepoint <- as.factor(dat_8_28_6800$Timepoint)
dat_8_28_6800$Ci <- as.numeric(dat_8_28_6800$Ci)

dat_8_28_600 <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/Li600_diurnal_august_filled.csv", header = FALSE)
dat_8_28_600 <- row_to_names(dat_8_28_600, 1)

dat_8_28_600 <- subset(dat_8_28_600, select = c(Time,Date,ID,gsw,VPDleaf,Fs,Fm,PhiPS2,ETR,Qamb,Tleaf))
str(dat_8_28_600)

#fix time format and add timepoint
?as.POSIXct
dat_8_28_600$Time <-format(strptime(dat_8_28_600$Time, "%I:%M:%S %p"), format="%H:%M:%S")
str(dat_8_28_600)

dat_8_28_600$Time <- hms(dat_8_28_600$Time)
dat_8_28_600$hour <- hour(dat_8_28_600$Time)
dat_8_28_600 <- dat_8_28_600 %>% mutate(Timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 16 & hour < 18 ~ 16,
  hour >= 18 & hour < 20 ~ 18,
  hour >= 20 & hour < 24 ~ 20
))

dat_8_28_600$gsw <- as.numeric(dat_8_28_600$gsw)
dat_8_28_600$VPDleaf <- as.numeric(dat_8_28_600$VPDleaf)
dat_8_28_600$Fs <- as.numeric(dat_8_28_600$Fs)
dat_8_28_600$Fm <- as.numeric(dat_8_28_600$Fm)
dat_8_28_600$PhiPS2 <- as.numeric(dat_8_28_600$PhiPS2)
dat_8_28_600$ETR <- as.numeric(dat_8_28_600$ETR)
dat_8_28_600$Qamb <- as.numeric(dat_8_28_600$Qamb)
dat_8_28_600$Tleaf <- as.numeric(dat_8_28_600$Tleaf)
dat_8_28_600$Timepoint <- as.factor(dat_8_28_600$Timepoint)

#dat_8_28 <- left_join(dat_8_28_600, dat_8_28_6800, by = c("ID","Timepoint"))

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))

dat_8_28_6800 <- inner_join(dat_8_28_6800, LC_meta, by = "ID")
dat_8_28_600 <- inner_join(dat_8_28_600, LC_meta, by = "ID")

dat_8_28_6800 <- dat_8_28_6800 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

dat_8_28_600 <- dat_8_28_600 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

#define control trees and drought trees 
control_trees <- read.csv(file="LC_2023/2023_Sampling/Other/Lc_irrigation_sampleII.csv")
drought_trees <- LC_meta[!(LC_meta$ID %in% control_trees$ID),]


dat_8_28_6800 <- dat_8_28_6800 %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))


dat_8_28_600 <- dat_8_28_600 %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))

str(dat_8_28_6800)
dat_8_28_6800$A <- as.numeric(dat_8_28_6800$A)
dat_8_28_6800$A_dark <- as.numeric(dat_8_28_6800$A_dark)

dat_8_28_6800$Qin <- as.numeric(dat_8_28_6800$Qin)
dat_8_28_6800$Tleaf <- as.numeric(dat_8_28_6800$Tleaf)

August_6800_summary_event <- subset(dat_8_28_6800, A > 0) %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
)

write.csv(August_6800_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_event.csv)")


August_6800_summary_tree <- subset(dat_8_28_6800, A > 0) %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  Fs = mean(Fs),
  Fm_prime = mean(Fm.),
  Fo. = mean(Fo.),
  A_dark = mean(A_dark),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
  
)

write.csv(August_6800_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_tree.csv)")

August_6800_summary_Class <- subset(dat_8_28_6800, A > 0) %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = (A_sd/(sqrt(n))),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
)

write.csv(August_6800_summary_Class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_class.csv)")


August_600_summary_event <- subset(dat_8_28_600,) %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n)))
)

write.csv(August_600_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_event.csv)")



August_600_summary_tree <- subset(dat_8_28_600,) %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n)))
)

write.csv(August_600_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_tree.csv)")

August_600_summary_class <- subset(dat_8_28_600,) %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n)))
)

write.csv(August_600_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_class.csv)")

#read in pre_dawn_data
#fluorescence 
Aug_pre_dawn_fluor <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/8_29_night_sampling_filled.csv")
Aug_pre_dawn_WP <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_field_conditions/August_diurnal_mPA.csv")

Aug_predawn_fluor_averaged <- Aug_pre_dawn_fluor %>% dplyr::group_by(ID) %>% dplyr::summarize(
  FvFm = mean(Fv.Fm),
  Fo = mean(Fo),
  Fm = mean(Fm)
)

Aug_predawn <- inner_join(Aug_predawn_fluor_averaged, Aug_pre_dawn_WP, by = "ID")

Aug_predawn <- inner_join(Aug_predawn, LC_meta, by = "ID")
Aug_predawn <- Aug_predawn %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))
Aug_predawn <- Aug_predawn %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))

Aug_predawn$Timepoint <- "6"

Aug_predawn_summary_event <- Aug_predawn %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  FvFm_sd = sd(FvFm),
  FvFm = mean(FvFm),
  FvFm_se = FvFm_sd/(sqrt(n)),
  Fo_sd = sd(Fo),
  Fo = mean(Fo),
  Fo_se = Fo_sd/(sqrt(n)),
  Fm_sd = sd(Fm),
  Fm = mean(Fm),
  Fm_se = Fm_sd/(sqrt(n)),
  mPA_sd = sd(mPA),
  mPA = mean(mPA),
  mPA_se = mPA_sd/(sqrt(n))
)

write.csv(Aug_predawn_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_event.csv)")

Aug_predawn_summary_tree <- Aug_predawn %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  FvFm_sd = sd(FvFm),
  FvFm = mean(FvFm),
  FvFm_se = FvFm_sd/(sqrt(n)),
  Fo_sd = sd(Fo),
  Fo = mean(Fo),
  Fo_se = Fo_sd/(sqrt(n)),
  Fm_sd = sd(Fm),
  Fm = mean(Fm),
  Fm_se = Fm_sd/(sqrt(n)),
  mPA_sd = sd(mPA),
  mPA = mean(mPA),
  mPA_se = mPA_sd/(sqrt(n))
)

write.csv(Aug_predawn_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_tree.csv)")

Aug_predawn_summary_class <- Aug_predawn %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
  n = n(),
  FvFm_sd = sd(FvFm),
  FvFm = mean(FvFm),
  FvFm_se = FvFm_sd/(sqrt(n)),
  Fo_sd = sd(Fo),
  Fo = mean(Fo),
  Fo_se = Fo_sd/(sqrt(n)),
  Fm_sd = sd(Fm),
  Fm = mean(Fm),
  Fm_se = Fm_sd/(sqrt(n)),
  mPA_sd = sd(mPA),
  mPA = mean(mPA),
  mPA_se = mPA_sd/(sqrt(n))
)

write.csv(Aug_predawn_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_class.csv)")


#Detailed photosynthetic parameters for subset of trees
head(August_6800_summary_tree)
head(Aug_predawn_summary_tree)
Aug_predawn_summary_tree$Timepoint <- as.factor(Aug_predawn_summary_tree$Timepoint)

August_fluor_set <- inner_join(August_6800_summary_tree, Aug_predawn_summary_tree, by = c("ID","treatment","n"))
August_fluor_set$Timepoint <- August_fluor_set$Timepoint.x
#introduce ETR and PAR from li600
August_600_subset <-subset(August_600_summary_tree, select = c(ID,Timepoint,treatment,ETR,PAR,Tleaf,PhiPS2,gsw))

August_fluor_set <- inner_join(August_fluor_set, August_600_subset, by = c("ID","treatment","Timepoint") )

#introduce respiration rates
August_dark_resp <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_filled/August_resp.csv")
August_fluor_set <- left_join(August_fluor_set, August_dark_resp)

#impute missing Rd measurements from assimilation
#the culprit
August_fluor_set[53:56,31]

resp_imput_set <- August_fluor_set %>% group_by(ID) %>% summarize(
  A = mean(A),
  Rd = mean(Rd)
)

Resp_imput <- lm(Rd~A, data=resp_imput_set[-12,])
plot(Resp_imput, which = 1)
plot(Resp_imput, which = 2)
plot(Resp_imput, which = 3)
plot(Resp_imput, which = 4)
plot(Resp_imput, which = 5)
#looks "ok" to me
ggplot(resp_imput_set,aes(x=A,y=Rd))+
  geom_point()
# a strong outlier, but there is a relationship. As mean assimilation increases, respiration rate becomes more negative as expected by theory.
summary(Resp_imput)
#rd = -0.06084*Amean - 0.62210
predict(Resp_imput, resp_imput_set, interval = "prediction")


August_fluor_set[53:56,43] <- -1.0730629

#impute temperature at dark respiration measurement from gsw

resptemp_imput_set <- August_fluor_set %>% group_by(ID) %>% summarize(
  gsw = mean(gsw),
  Tleaf = mean(Tleaf)
)

Resptemp_imput <- lm(Tleaf~gsw, data=resptemp_imput_set[-12,])
plot(Resptemp_imput, which = 1)
plot(Resptemp_imput, which = 2)
plot(Resptemp_imput, which = 3)
plot(Resptemp_imput, which = 4)
plot(Resptemp_imput, which = 5)
#looks "ok" to me
ggplot(resptemp_imput_set,aes(x=gsw,y=Tleaf))+
  geom_point()
#but there is a relationship
summary(Resptemp_imput)
#rd = -0.06084*Amean - 0.62210
predict(Resptemp_imput, resptemp_imput_set, interval = "prediction")

August_fluor_set[53:56,46]
August_fluor_set[53:56,46] <- 21.35103

#parameters derived from li600 for ETR, PAR and Tleaf
August_fluor_set$R <- (-1*August_fluor_set$Rd*2.2)^((August_fluor_set$Tleaf - August_fluor_set$Resp_temp)/10)

(August_fluor_set$Rd*2.2)^August_fluor_set$resp_factor
August_fluor_set$resp_factor <- (August_fluor_set$Tleaf - August_fluor_set$Resp_temp)/10

August_fluor_set$Theta_e <- (4*(August_fluor_set$PhiPS2 - 0.026947))/7.567127
August_fluor_set$Jt <- August_fluor_set$Theta_e*August_fluor_set$PAR
August_fluor_set$Jo <- (2/3)*(August_fluor_set$Jt - (4*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Rl <- (1/12)*(August_fluor_set$Jt - (4*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Jc <- (1/3)*(August_fluor_set$Jt + (8*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Jo_ratio <- August_fluor_set$Jo/August_fluor_set$Jt

#parameters for ETR, phips2 and Qin from 6800, tleaf from 600

August_fluor_set$R_6800 <- (-1*August_fluor_set$Rd*2.2)^((August_fluor_set$Tleaf - August_fluor_set$Resp_temp)/10)

August_fluor_set$Theta_e_6800 <- (4*(August_fluor_set$PhiPS2_6800 - 0.026947))/7.567127
#values derived from low O2 conditions
August_fluor_set$Jt_6800 <- August_fluor_set$Theta_e_6800*August_fluor_set$Qin
August_fluor_set$Jo_6800 <- (2/3)*(August_fluor_set$Jt_6800 - (4*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Rl_6800 <- (1/12)*(August_fluor_set$Jt_6800 - (4*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Jc_6800 <- (1/3)*(August_fluor_set$Jt_6800 + (8*(August_fluor_set$A + August_fluor_set$Rd)))
August_fluor_set$Jo_ratio_6800 <- August_fluor_set$Jo_6800/August_fluor_set$Jt_6800

August_fluor_set <- inner_join(August_fluor_set, LC_meta, by = "ID")

August_fluor_set <- August_fluor_set %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

August_fluor_set <- August_fluor_set %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))





August_fluor_set[,43]

August_fluor_set_clean <- subset(August_fluor_set, !is.na(August_fluor_set[,43]))

August_fluor_set_summary_event <- August_fluor_set_clean %>% dplyr::group_by(event_short,Timepoint,treatment) %>% dplyr::summarize(
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
  gsw_se = (gsw_sd/(sqrt(n))),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Jt_6800_sd = sd(Jt_6800),
  Jt_6800 = mean(Jt_6800),
  Jt_6800_se = (Jt_6800_sd/(sqrt(n))),
  Jo_6800_sd = sd(Jo_6800),
  Jo_6800 = mean(Jo_6800),
  Jo_6800_se = (Jo_6800_sd/(sqrt(n))),
  Jc_6800_sd = sd(Jc_6800),
  Jc_6800 = mean(Jc_6800),
  Jc_6800_se = (Jc_6800_sd/(sqrt(n))),
  Jo_ratio_6800_sd = sd(Jo_ratio_6800),
  Jo_ratio_6800 = mean(Jo_ratio_6800),
  Jo_ratio_6800_se = (Jo_ratio_6800_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
  
)  
 


write.csv(August_fluor_set_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_event.csv)")

August_fluor_set_summary_class <- August_fluor_set_clean %>% dplyr::group_by(Class,Timepoint,treatment) %>% dplyr::summarize(
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
  gsw_se = (gsw_sd/(sqrt(n))),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Jt_6800_sd = sd(Jt_6800),
  Jt_6800 = mean(Jt_6800),
  Jt_6800_se = (Jt_6800_sd/(sqrt(n))),
  Jo_6800_sd = sd(Jo_6800),
  Jo_6800 = mean(Jo_6800),
  Jo_6800_se = (Jo_6800_sd/(sqrt(n))),
  Jc_6800_sd = sd(Jc_6800),
  Jc_6800 = mean(Jc_6800),
  Jc_6800_se = (Jc_6800_sd/(sqrt(n))),
  Jo_ratio_6800_sd = sd(Jo_ratio_6800),
  Jo_ratio_6800 = mean(Jo_ratio_6800),
  Jo_ratio_6800_se = (Jo_ratio_6800_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
)

write.csv(August_fluor_set_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_class.csv)")

August_fluor_set_summary_tree <- August_fluor_set_clean %>% dplyr::group_by(ID,Timepoint,treatment) %>% dplyr::summarize(
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
  gsw_se = (gsw_sd/(sqrt(n))),
  ETR_6800_sd = sd(ETR),
  ETR_6800 = mean(ETR),
  ETR_6800_se = (ETR_6800_sd/(sqrt(n))),
  Tleaf_6800_sd = sd(Tleaf),
  Tleaf_6800 = mean(Tleaf),
  Tleaf_6800_se = (Tleaf_6800_sd/(sqrt(n))),
  PhiPS2_6800_sd = sd(PhiPS2),
  PhiPS2_6800 = mean(PhiPS2),
  PhiPS2_6800_se = (PhiPS2_6800_sd/(sqrt(n))),
  Qin_sd = sd(Qin),
  Qin = mean(Qin),
  Qin_se = (Qin_sd/(sqrt(n))),
  Jt_6800_sd = sd(Jt_6800),
  Jt_6800 = mean(Jt_6800),
  Jt_6800_se = (Jt_6800_sd/(sqrt(n))),
  Jo_6800_sd = sd(Jo_6800),
  Jo_6800 = mean(Jo_6800),
  Jo_6800_se = (Jo_6800_sd/(sqrt(n))),
  Jc_6800_sd = sd(Jc_6800),
  Jc_6800 = mean(Jc_6800),
  Jc_6800_se = (Jc_6800_sd/(sqrt(n))),
  Jo_ratio_6800_sd = sd(Jo_ratio_6800),
  Jo_ratio_6800 = mean(Jo_ratio_6800),
  Jo_ratio_6800_se = (Jo_ratio_6800_sd/(sqrt(n))),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = (Ci_sd/(sqrt(n)))
)

write.csv(August_fluor_set_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_tree.csv)")


