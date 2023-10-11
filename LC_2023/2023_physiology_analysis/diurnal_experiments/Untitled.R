####Analysis of diurnal curves on 6/27######


library(treenetproc)
library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(tidyr)

dat_6_27_6800 <- read.csv(file = "LC_2023/LiCOR/diurnal_experiments/6_27_diurnal_datasheet_li6800.csv")

dat_6_27_6800 <- subset(dat_6_27_6800, select = c(date,ID,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf))
