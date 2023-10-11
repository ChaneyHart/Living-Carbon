##Metabolite phys

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


July_phys_6800_1 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/Li6800/July_12-metabolite_6800.csv")
July_phys_6800_2 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/Li6800/July_13_metabolite_6800.csv")

July_phys_6800 <- rbind(July_phys_6800_1, July_phys_6800_2)

#fix time stamp
stri_sub(July_phys_6800$date,1,8)
stri_sub(July_phys_6800$date,5,4) <- "/"
stri_sub(July_phys_6800$date,8,7) <- "/"
July_phys_6800$date
July_phys_6800$date <- gsub('\\/',"-",July_phys_6800$date)
July_phys_6800$date <- ymd_hms(July_phys_6800$date)

July_phys_6800 <- subset(July_phys_6800, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf,A_dark))
str(July_phys_6800)
July_phys_6800$PhiPS2 <- as.numeric(July_phys_6800$PhiPS2)
July_phys_6800$ETR <- as.numeric(July_phys_6800$ETR)
July_phys_6800$NPQ <- as.numeric(July_phys_6800$NPQ)
July_phys_6800$Fo. <- as.numeric(July_phys_6800$Fo.)
July_phys_6800$qP <- as.numeric(July_phys_6800$qP)
July_phys_6800$qN <- as.numeric(July_phys_6800$qN)
July_phys_6800$qL <- as.numeric(July_phys_6800$qL)

#li600

July_phys_600_1 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/Li600/2023-07-12/LICOR_7_12_23_metabolite_phys_1_600.csv")
July_phys_600_1 <- July_phys_600_1[-c(2),]
July_phys_600_1 <- row_to_names(July_phys_600_1, 1)
July_phys_600_2 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/Li600/2023-07-13/LICOR_7_13_23_metabolite_phys_2_600.csv")
July_phys_600_2 <- July_phys_600_2[-c(2),]
July_phys_600_2 <- row_to_names(July_phys_600_2, 1)

July_phys_600 <- rbind(July_phys_600_1, July_phys_600_2)
July_phys_600 <- subset(July_phys_600, select = c(Time,Date,ID,gsw,VPDleaf,Fs,Fm,PhiPS2,ETR,Qamb,Tleaf))
str(July_phys_600)

July_phys_600$gsw <- as.numeric(July_phys_600$gsw)
July_phys_600$VPDleaf <- as.numeric(July_phys_600$VPDleaf)
July_phys_600$Fs <- as.numeric(July_phys_600$Fs)
July_phys_600$Fm <- as.numeric(July_phys_600$Fm)
July_phys_600$PhiPS2 <- as.numeric(July_phys_600$PhiPS2)
July_phys_600$ETR <- as.numeric(July_phys_600$ETR)
July_phys_600$Qamb <- as.numeric(July_phys_600$Qamb)
July_phys_600$Tleaf <- as.numeric(July_phys_600$Tleaf)
July_phys_600$ETR <- as.numeric(July_phys_600$ETR)


July_phys <- inner_join(July_phys_6800, July_phys_600, by = "ID")

#August




