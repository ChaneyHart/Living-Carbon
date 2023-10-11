#Compile response curves and break into Light response and CO2 response datasets
#data collected summer 2023


library(plantecophys)
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(plantecophys)
library(mgcv)

#read in individual datasheets and compile
str(curves_6_14_2)
#June 14 & 15
curves_6_14_1 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_6_14_1_proc.csv")
curves_6_14_1 <- subset(curves_6_14_1, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_6_14_2 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_6_14_2_proc.csv")
curves_6_14_2 <- subset(curves_6_14_2, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_6_15 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_6_15_proc.csv")
curves_6_15 <- subset(curves_6_15, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_Jun14 <- rbind(curves_6_14_1, curves_6_14_2, curves_6_15)
curves_Jun14$sample_date <- "Jun14"

#June 22 & 23
curves_6_22 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_6_22_proc.csv")
curves_6_22 <- subset(curves_6_22, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_6_23 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_6_23_proc.csv")
curves_6_23 <- subset(curves_6_23, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_Jun22 <- rbind(curves_6_22, curves_6_23)
curves_Jun22$sample_date <- "Jun22"

#July 18 & 19

curves_7_18 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_18_proc.csv")
curves_7_18 <- subset(curves_7_18, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_7_19 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_19_proc.csv")
curves_7_19 <- subset(curves_7_19, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_July18 <- rbind(curves_7_18, curves_7_19)
curves_July18$sample_date <- "July18"

#July 20 & 21

curves_7_20 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_20_proc.csv")
curves_7_20 <- subset(curves_7_20, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_7_21 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_21_proc.csv")
curves_7_21 <- subset(curves_7_21, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_July20 <- rbind(curves_7_20, curves_7_21)
curves_July20$sample_date <- "July20"

#July 27 & 28

curves_7_27 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_27_proc.csv")
curves_7_27 <- subset(curves_7_27, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_7_28 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_7_28_proc.csv")
curves_7_28 <- subset(curves_7_28, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_July27 <- rbind(curves_7_27, curves_7_28)
curves_July27$sample_date <- "July27"


#August 25 and 29

curves_8_25 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_8_25_proc.csv")
curves_8_25 <- subset(curves_8_25, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_8_29_1 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_8_29_proc_1.csv")
curves_8_29_1 <- subset(curves_8_29_1, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_8_29_2 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_8_29_proc_2.csv")
curves_8_29_2 <- subset(curves_8_29_2, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_8_29_3 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/proc_csvs/dat_8_29_proc_3.csv")
curves_8_29_3 <- subset(curves_8_29_3, select = c(date,SPAD,Midday_WP,ID,Response_curve,E,Emm,A,Ca,Ci,CO2_r_sp,Qin,gsw,RHcham,VPcham,VPDleaf,Fm,Fv.Fm,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,Fv..Fm.,qP,qN,qP_Fo,qN_Fo,Tleaf))

curves_Aug29 <- rbind(curves_8_25, curves_8_29_1, curves_8_29_2, curves_8_29_3)
curves_Aug29$sample_date <- "Aug29"

####Compile datasheets from seperate dates together

response_curves_compiled <- rbind(curves_Jun14,curves_Jun22,curves_July18,curves_July20,curves_July27,curves_Aug29)

light_response_compiled <- subset(response_curves_compiled, Response_curve == "Q")
write.csv(light_response_compiled, file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/light_response_curves_compiled.csv")
CO2_response_compiled <- subset(response_curves_compiled, Response_curve == "CO2")
write.csv(CO2_response_compiled, file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/CO2_response_curves_compiled.csv")

