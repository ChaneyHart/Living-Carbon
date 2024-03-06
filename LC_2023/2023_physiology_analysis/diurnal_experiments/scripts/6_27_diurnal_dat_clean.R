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
dat_6_27_6800$Ci <- as.numeric(dat_6_27_6800$Ci)
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
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

dat_6_27_600 <- dat_6_27_600 %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

#Create Summary table and graph


June_6800_summary_event <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(event_short,Timepoint) %>% dplyr::summarize(
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

write.csv(June_6800_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_event.csv)")

June_6800_summary_tree <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(ID,Timepoint) %>% dplyr::summarize(
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n)))
  
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n))),
  
)

write.csv(June_600_summary_tree, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_tree.csv)")


June_6800_summary_class <- subset(dat_6_27_6800, A > 0) %>% dplyr::group_by(Class,Timepoint) %>% dplyr::summarize(
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
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = (PhiPS2_sd/(sqrt(n)))
)

write.csv(June_600_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_class.csv)")

# join datasets together

#June_summary <- inner_join(June_600_summary, June_6800_summary, by = c("event_short","Timepoint"))

June_summary_class <- inner_join(June_600_summary_class, June_6800_summary_class, by = c("Class","Timepoint"))

#respiration data

June_resp_dat <- subset(dat_6_27_6800, A < 0)
June_resp_dat$Rd <- June_resp_dat$A
June_resp_dat$Resp_temp <- June_resp_dat$Tleaf

June_resp_dat <- subset(June_resp_dat, select = c(ID,Rd,Resp_temp))
#join with data that has full measurements from 6800 and 600

June_fluor_set <- inner_join(June_6800_summary_tree,June_600_summary_tree, by = c("ID","Timepoint"))

June_fluor_set <- left_join(June_fluor_set, June_resp_dat)

June_fluor_set$R <- (-1*June_fluor_set$Rd*2.2)^((June_fluor_set$Tleaf - June_fluor_set$Resp_temp)/10)

June_fluor_set$Theta_e <- (4*(June_fluor_set$PhiPS2 - 0.026947))/7.567127
#values derived from low O2 conditions
June_fluor_set$Jt <- June_fluor_set$Theta_e*June_fluor_set$PAR
June_fluor_set$Jo <- (2/3)*(June_fluor_set$Jt - (4*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Rl <- (1/12)*(June_fluor_set$Jt - (4*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Jc <- (1/3)*(June_fluor_set$Jt + (8*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Jo_ratio <- June_fluor_set$Jo/June_fluor_set$Jt

#parameters derived from 6800 data for ETR and Tleaf
June_fluor_set$R_6800 <- (-1*June_fluor_set$Rd*2.2)^((June_fluor_set$Tleaf - June_fluor_set$Resp_temp)/10)

June_fluor_set$Theta_e_6800 <- (4*(June_fluor_set$PhiPS2_6800 - 0.026947))/7.567127
#values derived from low O2 conditions
June_fluor_set$Jt_6800 <- June_fluor_set$Theta_e_6800*June_fluor_set$Qin
June_fluor_set$Jo_6800 <- (2/3)*(June_fluor_set$Jt_6800 - (4*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Rl_6800 <- (1/12)*(June_fluor_set$Jt_6800 - (4*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Jc_6800 <- (1/3)*(June_fluor_set$Jt_6800 + (8*(June_fluor_set$A + June_fluor_set$Rd)))
June_fluor_set$Jo_ratio_6800 <- June_fluor_set$Jo_6800/June_fluor_set$Jt_6800



June_fluor_set <- inner_join(June_fluor_set, LC_meta, by = "ID")

June_fluor_set <- June_fluor_set %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT"))

June_fluor_set <- June_fluor_set %>% mutate(treatment = case_when(
  ID %in% control_trees$ID ~ "control",
  ID %in% drought_trees$ID ~ "drought"
))

June_fluor_set_clean <- filter(June_fluor_set, R > 0)
write.csv(June_fluor_set_clean, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_fluor_set_clean.csv")

June_fluor_set_summary_event <- June_fluor_set_clean %>% dplyr::group_by(event_short,Timepoint) %>% dplyr::summarize(
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



write.csv(June_fluor_set_summary_event, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_28_fluor_set_diurnal_summary_event.csv)")

June_fluor_set_summary_class <- June_fluor_set_clean %>% dplyr::group_by(Class,Timepoint) %>% dplyr::summarize(
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
  Rl_sd = sd(Rl),
  Rl = mean(Rl),
  Rl_se = (Rl_sd/(sqrt(n))),
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

write.csv(June_fluor_set_summary_class, file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_28_fluor_set_diurnal_summary_class.csv)")

June_fluor_set_summary_tree <- June_fluor_set_clean %>% dplyr::group_by(ID,Timepoint) %>% dplyr::summarize(
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


#graph
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

June_Jo_6800 <- ggplot(June_fluor_set_summary_class, aes(x = Timepoint, group=Class,color = Class))+
  geom_point(aes(y = Jo_6800), size = 2)+
  geom_line(aes(y = Jo_6800))+ 
  geom_errorbar(aes(ymin = Jo_6800-Jo_6800_se, ymax = Jo_6800+Jo_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")

June_Jo_6800

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
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio +Jo_ratio_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
June_Jo_ratio

June_Jo_ratio_6800 <- ggplot(subset(June_fluor_set_summary_class, Timepoint != 18), aes(x = Timepoint, group = Class,color = Class))+
  geom_point(aes(y = Jo_ratio_6800), size = 2)+
  geom_line(aes(y = Jo_ratio_6800))+ 
  geom_errorbar(aes(ymin = Jo_ratio_6800-Jo_ratio_6800_se, ymax = Jo_ratio_6800+Jo_ratio_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")+
  theme_bw()
June_Jo_ratio_6800


June_photorespiration_plot <- grid.arrange(June_Jo, June_Rl, June_Jo_ratio, nrow=3)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_photorespiration_plot.png",plot = June_photorespiration_plot, width = 8, height = 6, dpi = 300)

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

June_SM <- read.csv("LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_field_conditions/6_27_moisture_data (filled).csv")
mean(June_SM$SM.Av)
mean(June_SM$SM.Av) - 2*(sd(June_SM$SM.Av)/(sqrt(48)))
mean(June_SM$SM.Av) + 2*(sd(June_SM$SM.Av)/(sqrt(48)))


#Bring in weather data
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
