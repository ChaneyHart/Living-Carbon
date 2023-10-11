####Analysis of diurnal curves on 6/16######


library(treenetproc)
library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(tidyr)

dat_6_15 <- read.csv(file = "LC_2023/LiCOR/diurnal_experiments/dat_6_16_diurnal.csv")

diurnal_6_15 <- subset(dat_6_15, select = c(date,tree,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf))

#Read in respiration dat
resp_dat_6_15 <- read.csv("LC_2023/LiCOR/diurnal_experiments/dat_6_15_respiration.csv")
resp_6_15 <- subset(resp_dat_6_15, select = c(date,tree,E,A,Ca,Ci,gsw,gbw,gtw,VPDleaf,Fs,Fm.,PhiPS2,ETR,PhiCO2,NPQ,Fmin,Fo.,qP,qN,qP_Fo,qN_Fo,qL,Qin,Tleaf))
#join
diurnal_6_15 <- rbind(diurnal_6_15, resp_6_15)

#fix ID label
diurnal_6_15$tree_2 <- substr(diurnal_6_15$tree,5,8)
diurnal_6_15$tree_1 <- toupper(substr(diurnal_6_15$tree,1,4))
diurnal_6_15$ID <- paste(diurnal_6_15$tree_1,diurnal_6_15$tree_2, sep = "")
#fix time stamp
diurnal_6_15$date <- gsub('\\/',"-",diurnal_6_15$date)
diurnal_6_15$date <- mdy_hm(diurnal_6_15$date)


#add meta data
LC_meta <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_inventory, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))

diurnal_6_15 <- inner_join(diurnal_6_15, LC_meta, by = "ID")

str(diurnal_6_15$date)

#add info about time of sampling
diurnal_6_15$hour <- hour(diurnal_6_15$date)


diurnal_6_15 <- diurnal_6_15 %>% mutate(timepoint = case_when(
  hour >= 8 & hour < 10 ~ 8,
  hour >= 10 & hour < 12 ~ 10,
  hour >= 12 & hour < 14 ~ 12,
  hour >= 14 & hour < 16 ~ 14,
  hour >= 15 & hour < 18 ~ 16,
  hour >= 20 & hour < 24 ~ 20
))

diurnal_6_15$timepoint <- as.factor(diurnal_6_15$timepoint)

##diurnal curves

ggplot(diurnal_6_15,aes(x= timepoint, y=A, fill = construct2))+
  geom_boxplot()

ggplot(diurnal_6_15,aes(x= timepoint, y=gsw, fill = construct2))+
  geom_boxplot()

ggplot(diurnal_6_15,aes(x= timepoint, y=PhiPS2, fill = construct2))+
  geom_boxplot()

#8-10
ggplot(subset(diurnal_6_15, timepoint == 8), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 8), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()

#10
ggplot(subset(diurnal_6_15, timepoint == 10), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 10), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()

#12
ggplot(subset(diurnal_6_15, timepoint == 12), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 12), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()

#14
ggplot(subset(diurnal_6_15, timepoint == 14), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 14), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()

#16
ggplot(subset(diurnal_6_15, timepoint == 16), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 16), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()

#20
ggplot(subset(diurnal_6_15, timepoint == 20), aes(x=event_short, y=A, fill = construct2))+
  geom_boxplot()

ggplot(subset(diurnal_6_15, timepoint == 20), aes(x=construct2, y=A, fill = construct2))+
  geom_boxplot()



ggplot(diurnal_6_15,aes(x= timepoint, y=A, fill = construct2))+
  geom_boxplot()
