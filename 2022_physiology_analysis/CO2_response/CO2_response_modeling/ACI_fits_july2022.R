
#A_Ci_Curves
#data collected June and July of 2022

#adapted from Adam Sibley's script

library(plantecophys)
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(plantecophys)
library(mgcv)

#set working directory

#read in data
LC_spreadz = '2022_physiology_analysis/CO2_response/compiled/ACI_individual_tree.xlsx' 


#Function for cleaning up data for fits
#viewfit = function(fit){plot(y = fit$df$Ameas, x = fit$df$Ci);fit}
viewdat = function(dat){plot(y = dat$A, x = dat$Ci);abline(h = seq(0,20,1));dat}

#fx for plotting fits

plotfits = function(fit){
  jpeg(filename = paste('',deparse(substitute(fit)),'.jpeg',sep=''),width = 10,height=7,units='in',res = 100)
  plot(fit,main =deparse(substitute(fit)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel),labels = paste('Vcmax =',round(fit$pars[1,1],2)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel)-1,labels = paste('Jmax =',round(fit$pars[2,1],2)))
  dev.off()
}

#creating data frame to store data

#create separate dataframe for observations with fluoresence data 
#subset data for which fluorescence was on  
all_trees_df_flr <- subset(all_trees_df, Fm. > 100)

#convert to wide format
all_trees_wide <- pivot_wider(all_trees_df,id_cols = "tree", names_from = Csetpoint, values_from = c(A,gsw,Ci))
all_trees_flr_wide <- pivot_wider(all_trees_df_flr,id_cols = "tree", names_from = Csetpoint, values_from = c(A,gsw,Ci,PhiPS2,ETR,PhiCO2))



tree_list = list()
Vcmax_list = list()
Jmax_list = list()
Rd_list = list()


##look at fits, model and graph them and get outputs and assign to storage dataframe

LCOR578 = read.xlsx(LC_spreadz,'LCOR-578')
viewdat(LCOR578)
LCOR578 <- LCOR578[c(-1,-2,-7),]
LCOR578_fit = fitaci(LCOR578,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR578_fit)
summary(LCOR578_fit)
coef578 <- coef(LCOR578_fit)
LCOR578_fit$Photosyn(Ci=200)
coef578
tree_list <- append(tree_list, 'LCOR-578')
Vcmax_list <- append(Vcmax_list, coef578[1])
Jmax_list <- append(Jmax_list, coef578[2])
Rd_list <- append(Rd_list, coef578[3])



LCOR269 = read.xlsx(LC_spreadz,'LCOR-269')
viewdat(LCOR269)
LCOR269 <- LCOR269[c(-1,-2,-7),]
LCOR269_fit = fitaci(LCOR269,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR269_fit)
summary(LCOR269_fit)
coef269 <- coef(LCOR269_fit)
LCOR269_fit$Photosyn(Ci=200)
coef269
tree_list <- append(tree_list, 'LCOR-269')
Vcmax_list <- append(Vcmax_list, coef269[1])
Jmax_list <- append(Jmax_list, coef269[2])
Rd_list <- append(Rd_list, coef269[3])


LCOR211 = read.xlsx(LC_spreadz,'LCOR-211')
viewdat(LCOR211)
LCOR211 <- LCOR211[c(-1,-2,-7),]
LCOR211_fit = fitaci(LCOR211,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR211_fit)
summary(LCOR211_fit)
coef211 <- coef(LCOR211_fit)
LCOR211_fit$Photosyn(Ci=200)
coef211
tree_list <- append(tree_list, 'LCOR-211')
Vcmax_list <- append(Vcmax_list, coef211[1])
Jmax_list <- append(Jmax_list, coef211[2])
Rd_list <- append(Rd_list, coef211[3])



LCOR152 = read.xlsx(LC_spreadz,'LCOR-152')
viewdat(LCOR152)
LCOR152 <- LCOR152[c(-1,-2,-7),]
LCOR152_fit = fitaci(LCOR152,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR152_fit)
summary(LCOR152_fit)
coef152 <- coef(LCOR152_fit)
LCOR152_fit$Photosyn(Ci=200)
coef152
tree_list <- append(tree_list, 'LCOR-152')
Vcmax_list <- append(Vcmax_list, coef152[1])
Jmax_list <- append(Jmax_list, coef152[2])
Rd_list <- append(Rd_list, coef152[3])

LCOR566 = read.xlsx(LC_spreadz,'LCOR-566')
viewdat(LCOR566)
LCOR566 <- LCOR566[c(-1,-2,-8),]
LCOR566_fit = fitaci(LCOR566,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR566_fit)
summary(LCOR566_fit)
coef566 <- coef(LCOR566_fit)
coef566
tree_list <- append(tree_list, 'LCOR-566')
Vcmax_list <- append(Vcmax_list, coef566[1])
Jmax_list <- append(Jmax_list, coef566[2])
Rd_list <- append(Rd_list, coef566[3])



LCOR079 = read.xlsx(LC_spreadz,'LCOR-079')
viewdat(LCOR079)
LCOR079 <- LCOR079[c(-1,-2,-8),]
LCOR079_fit = fitaci(LCOR079,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR079_fit)
summary(LCOR079_fit)
coef079 <- coef(LCOR079_fit)
coef079
tree_list <- append(tree_list, 'LCOR-079')
Vcmax_list <- append(Vcmax_list, coef079[1])
Jmax_list <- append(Jmax_list, coef079[2])
Rd_list <- append(Rd_list, coef079[3])


LCOR318 = read.xlsx(LC_spreadz,'LCOR-318')
viewdat(LCOR318)
LCOR318 <- LCOR318[c(-1,-2,-8),]
LCOR318_fit = fitaci(LCOR318,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR318_fit)
summary(LCOR318_fit)
coef318 <- coef(LCOR318_fit)
tree_list <- append(tree_list, 'LCOR-318')
Vcmax_list <- append(Vcmax_list, coef318[1])
Jmax_list <- append(Jmax_list, coef318[2])
Rd_list <- append(Rd_list, coef318[3])


LCOR160 = read.xlsx(LC_spreadz,'LCOR-160')
viewdat(LCOR160)
LCOR160 <- LCOR160[c(-1,-2,-8),]
LCOR160_fit = fitaci(LCOR160,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR160_fit)
summary(LCOR160_fit)
coef160 <- coef(LCOR160_fit)
tree_list <- append(tree_list, 'LCOR-160')
Vcmax_list <- append(Vcmax_list, coef160[1])
Jmax_list <- append(Jmax_list, coef160[2])
Rd_list <- append(Rd_list, coef160[3])


LCOR106 = read.xlsx(LC_spreadz,'LCOR-106')
viewdat(LCOR106)
LCOR106 <- LCOR106[c(-1,-2,-8),]
LCOR106_fit = fitaci(LCOR106,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR106_fit)
summary(LCOR106_fit)
coef106 <- coef(LCOR106_fit)
tree_list <- append(tree_list, 'LCOR-106')
Vcmax_list <- append(Vcmax_list, coef106[1])
Jmax_list <- append(Jmax_list, coef106[2])
Rd_list <- append(Rd_list, coef106[3])


LCOR509 = read.xlsx(LC_spreadz,'LCOR-509')
viewdat(LCOR509)
LCOR509 <- LCOR509[c(-1,-2,-8),]
LCOR509_fit = fitaci(LCOR509,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR509_fit)
summary(LCOR509_fit)
coef509 <- coef(LCOR509_fit)
tree_list <- append(tree_list, 'LCOR-509')
Vcmax_list <- append(Vcmax_list, coef509[1])
Jmax_list <- append(Jmax_list, coef509[2])
Rd_list <- append(Rd_list, coef509[3])


LCOR474 = read.xlsx(LC_spreadz,'LCOR-474')
viewdat(LCOR474)
LCOR474 <- LCOR474[c(-1,-2,-8),]
LCOR474_fit = fitaci(LCOR474,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR474_fit)
summary(LCOR474_fit)
coef474 <- coef(LCOR474_fit)
tree_list <- append(tree_list, 'LCOR-474')
Vcmax_list <- append(Vcmax_list, coef474[1])
Jmax_list <- append(Jmax_list, coef474[2])
Rd_list <- append(Rd_list, coef474[3])


LCOR217 = read.xlsx(LC_spreadz,'LCOR-217')
viewdat(LCOR217)
LCOR217 <- LCOR217[c(-1,-2,-8),]
LCOR217_fit = fitaci(LCOR217,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR217_fit)
summary(LCOR217_fit)
coef217 <- coef(LCOR217_fit)
tree_list <- append(tree_list, 'LCOR-217')
Vcmax_list <- append(Vcmax_list, coef217[1])
Jmax_list <- append(Jmax_list, coef217[2])
Rd_list <- append(Rd_list, coef217[3])


LCOR231 = read.xlsx(LC_spreadz,'LCOR-231')
viewdat(LCOR231)
LCOR231 <- LCOR231[c(-1,-2,-8),]
LCOR231_fit = fitaci(LCOR231,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR231_fit)
summary(LCOR231_fit)
coef231 <- coef(LCOR231_fit)
tree_list <- append(tree_list, 'LCOR-231')
Vcmax_list <- append(Vcmax_list, coef231[1])
Jmax_list <- append(Jmax_list, coef231[2])
Rd_list <- append(Rd_list, coef231[3])


LCOR108 = read.xlsx(LC_spreadz,'LCOR-108')
viewdat(LCOR108)
LCOR108 <- LCOR108[c(-1,-2,-8),]
LCOR108_fit = fitaci(LCOR108,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR108_fit)
summary(LCOR108_fit)
coef108 <- coef(LCOR108_fit)
tree_list <- append(tree_list, 'LCOR-108')
Vcmax_list <- append(Vcmax_list, coef108[1])
Jmax_list <- append(Jmax_list, coef108[2])
Rd_list <- append(Rd_list, coef108[3])


LCOR584 = read.xlsx(LC_spreadz,'LCOR-584_usetleafeb')
viewdat(LCOR584)
LCOR584 <- LCOR584[c(-1,-2,-8),]
LCOR584_fit = fitaci(LCOR584,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR584_fit)
summary(LCOR584_fit)
coef584 <- coef(LCOR584_fit)
tree_list <- append(tree_list, 'LCOR-584')
Vcmax_list <- append(Vcmax_list, coef584[1])
Jmax_list <- append(Jmax_list, coef584[2])
Rd_list <- append(Rd_list, coef584[3])


LCOR225 = read.xlsx(LC_spreadz,'LCOR-225_usetleafeb')
viewdat(LCOR225)
LCOR225 <- LCOR225[c(-1,-2,-8),]
LCOR225_fit = fitaci(LCOR225,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR225_fit)
summary(LCOR225_fit)
coef225 <- coef(LCOR225_fit)
tree_list <- append(tree_list, 'LCOR-225')
Vcmax_list <- append(Vcmax_list, coef225[1])
Jmax_list <- append(Jmax_list, coef225[2])
Rd_list <- append(Rd_list, coef225[3])


LCOR472 = read.xlsx(LC_spreadz,'LCOR-472_usetleafeb')
viewdat(LCOR472)
LCOR472 <- LCOR472[c(-1,-2,-8),]
LCOR472_fit = fitaci(LCOR472,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR472_fit)
summary(LCOR472_fit)
coef472 <- coef(LCOR472_fit)
tree_list <- append(tree_list, 'LCOR-472')
Vcmax_list <- append(Vcmax_list, coef472[1])
Jmax_list <- append(Jmax_list, coef472[2])
Rd_list <- append(Rd_list, coef472[3])


LCOR586 = read.xlsx(LC_spreadz,'LCOR-586_usetleafeb')
viewdat(LCOR586)
LCOR586 <- LCOR586[c(-1,-2,-8),]
LCOR586_fit = fitaci(LCOR586,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR586_fit)
summary(LCOR586_fit)
coef586 <- coef(LCOR586_fit)
tree_list <- append(tree_list, 'LCOR-586')
Vcmax_list <- append(Vcmax_list, coef586[1])
Jmax_list <- append(Jmax_list, coef586[2])
Rd_list <- append(Rd_list, coef586[3])


LCOR383 = read.xlsx(LC_spreadz,'LCOR-383')
viewdat(LCOR383)
LCOR383 <- LCOR383[c(-1,-2,-8),]
LCOR383_fit = fitaci(LCOR383,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR383_fit)
summary(LCOR383_fit)
coef383 <- coef(LCOR383_fit)
tree_list <- append(tree_list, 'LCOR-383')
Vcmax_list <- append(Vcmax_list, coef383[1])
Jmax_list <- append(Jmax_list, coef383[2])
Rd_list <- append(Rd_list, coef383[3])


LCOR001 = read.xlsx(LC_spreadz,'LCOR-001')
viewdat(LCOR001)
LCOR001 <- LCOR001[c(-1,-2,-8),]
LCOR001_fit = fitaci(LCOR001,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR001_fit)
summary(LCOR001_fit)
coef001 <- coef(LCOR001_fit)
tree_list <- append(tree_list, 'LCOR-001')
Vcmax_list <- append(Vcmax_list, coef001[1])
Jmax_list <- append(Jmax_list, coef001[2])
Rd_list <- append(Rd_list, coef001[3])


LCOR546 = read.xlsx(LC_spreadz,'LCOR-546')
viewdat(LCOR546)
LCOR546 <- LCOR546[c(-1,-2,-8),]
LCOR546_fit = fitaci(LCOR546,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR546_fit)
summary(LCOR546_fit)
coef546 <- coef(LCOR546_fit)
tree_list <- append(tree_list, 'LCOR-546')
Vcmax_list <- append(Vcmax_list, coef546[1])
Jmax_list <- append(Jmax_list, coef546[2])
Rd_list <- append(Rd_list, coef546[3])


LCOR501 = read.xlsx(LC_spreadz,'LCOR-501')
viewdat(LCOR501)
LCOR501 <- LCOR501[c(-1,-2,-8),]
LCOR501_fit = fitaci(LCOR501,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR501_fit)
summary(LCOR501_fit)
coef501 <- coef(LCOR501_fit)
tree_list <- append(tree_list, 'LCOR-501')
Vcmax_list <- append(Vcmax_list, coef501[1])
Jmax_list <- append(Jmax_list, coef501[2])
Rd_list <- append(Rd_list, coef501[3])


LCOR219 = read.xlsx(LC_spreadz,'LCOR-219_usetleafeb')
viewdat(LCOR219)
LCOR219 <- LCOR219[c(-1,-2,-8),]
LCOR219_fit = fitaci(LCOR219,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR219_fit)
summary(LCOR219_fit)
coef219 <- coef(LCOR219_fit)
tree_list <- append(tree_list, 'LCOR-219')
Vcmax_list <- append(Vcmax_list, coef219[1])
Jmax_list <- append(Jmax_list, coef219[2])
Rd_list <- append(Rd_list, coef219[3])


LCOR033 = read.xlsx(LC_spreadz,'LCOR-033_usetleafeb')
viewdat(LCOR033)
LCOR033 <- LCOR033[c(-1,-2,-8),]
LCOR033_fit = fitaci(LCOR033,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'TleafEB',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR033_fit)
summary(LCOR033_fit)
coef033 <- coef(LCOR033_fit)
tree_list <- append(tree_list, 'LCOR-033')
Vcmax_list <- append(Vcmax_list, coef033[1])
Jmax_list <- append(Jmax_list, coef033[2])
Rd_list <- append(Rd_list, coef033[3])


LCOR187 = read.xlsx(LC_spreadz,'LCOR-187')
viewdat(LCOR187)
LCOR187 <- LCOR187[c(-1,-2,-8),]
LCOR187_fit = fitaci(LCOR187,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR187_fit)
summary(LCOR187_fit)
coef187 <- coef(LCOR187_fit)
tree_list <- append(tree_list, 'LCOR-187')
Vcmax_list <- append(Vcmax_list, coef187[1])
Jmax_list <- append(Jmax_list, coef187[2])
Rd_list <- append(Rd_list, coef187[3])


LCOR243 = read.xlsx(LC_spreadz,'LCOR-243')
viewdat(LCOR243)
LCOR243 <- LCOR243[c(-1,-2,-8),]
LCOR243_fit = fitaci(LCOR243,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR243_fit)
summary(LCOR243_fit)
coef243 <- coef(LCOR243_fit)
tree_list <- append(tree_list, 'LCOR-243')
Vcmax_list <- append(Vcmax_list, coef243[1])
Jmax_list <- append(Jmax_list, coef243[2])
Rd_list <- append(Rd_list, coef243[3])


LCOR389 = read.xlsx(LC_spreadz,'LCOR-389')
viewdat(LCOR389)
LCOR389 <- LCOR389[c(-1,-2,-8),]
LCOR389_fit = fitaci(LCOR389,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR389_fit)
summary(LCOR389_fit)
coef389 <- coef(LCOR389_fit)
tree_list <- append(tree_list, 'LCOR-389')
Vcmax_list <- append(Vcmax_list, coef389[1])
Jmax_list <- append(Jmax_list, coef389[2])
Rd_list <- append(Rd_list, coef389[3])


LCOR392 = read.xlsx(LC_spreadz,'LCOR-392')
viewdat(LCOR392)
LCOR392 <- LCOR392[c(-1,-7),]
LCOR392_fit = fitaci(LCOR392,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR392_fit)
summary(LCOR392_fit)
coef392 <- coef(LCOR392_fit)
tree_list <- append(tree_list, 'LCOR-392')
Vcmax_list <- append(Vcmax_list, coef392[1])
Jmax_list <- append(Jmax_list, coef392[2])
Rd_list <- append(Rd_list, coef392[3])


LCOR500 = read.xlsx(LC_spreadz,'LCOR-500')
viewdat(LCOR500)
LCOR500 <- LCOR500[c(-1,-2,-8),]
LCOR500_fit = fitaci(LCOR500,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR500_fit)
summary(LCOR500_fit)
coef500 <- coef(LCOR500_fit)
tree_list <- append(tree_list, 'LCOR-500')
Vcmax_list <- append(Vcmax_list, coef500[1])
Jmax_list <- append(Jmax_list, coef500[2])
Rd_list <- append(Rd_list, coef500[3])


LCOR190 = read.xlsx(LC_spreadz,'LCOR-190')
viewdat(LCOR190)
LCOR190 <- LCOR190[c(-1,-2,-8),]
LCOR190_fit = fitaci(LCOR190,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR190_fit)
summary(LCOR190_fit)
coef190 <- coef(LCOR190_fit)
tree_list <- append(tree_list, 'LCOR-190')
Vcmax_list <- append(Vcmax_list, coef190[1])
Jmax_list <- append(Jmax_list, coef190[2])
Rd_list <- append(Rd_list, coef190[3])


LCOR077 = read.xlsx(LC_spreadz,'LCOR-077')
viewdat(LCOR077)
LCOR077 <- LCOR077[c(-1,-2,-8),]
LCOR077_fit = fitaci(LCOR077,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR077_fit)
summary(LCOR077_fit)
coef077 <- coef(LCOR077_fit)
tree_list <- append(tree_list, 'LCOR-077')
Vcmax_list <- append(Vcmax_list, coef077[1])
Jmax_list <- append(Jmax_list, coef077[2])
Rd_list <- append(Rd_list, coef077[3])


LCOR065 = read.xlsx(LC_spreadz,'LCOR-065')
viewdat(LCOR065)
LCOR065 <- LCOR065[c(-1,-2,-8),]
LCOR065_fit = fitaci(LCOR065,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR065_fit)
summary(LCOR065_fit)
coef065 <- coef(LCOR065_fit)
tree_list <- append(tree_list, 'LCOR-065')
Vcmax_list <- append(Vcmax_list, coef065[1])
Jmax_list <- append(Jmax_list, coef065[2])
Rd_list <- append(Rd_list, coef065[3])

#bit wonky but ok

LCOR292 = read.xlsx(LC_spreadz,'LCOR-292')
viewdat(LCOR292)
LCOR292 <- LCOR292[c(-1,-2,-8),]
LCOR292_fit = fitaci(LCOR292,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR292_fit)
summary(LCOR292_fit)
coef292 <- coef(LCOR292_fit)
tree_list <- append(tree_list, 'LCOR-292')
Vcmax_list <- append(Vcmax_list, coef292[1])
Jmax_list <- append(Jmax_list, coef292[2])
Rd_list <- append(Rd_list, coef292[3])


LCOR212 = read.xlsx(LC_spreadz,'LCOR-212')
viewdat(LCOR212)
LCOR212 <- LCOR212[c(-1,-2,-8),]
LCOR212_fit = fitaci(LCOR212,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR212_fit)
summary(LCOR212_fit)
coef212 <- coef(LCOR212_fit)
tree_list <- append(tree_list, 'LCOR-212')
Vcmax_list <- append(Vcmax_list, coef212[1])
Jmax_list <- append(Jmax_list, coef212[2])
Rd_list <- append(Rd_list, coef212[3])


LCOR010 = read.xlsx(LC_spreadz,'LCOR-010')
viewdat(LCOR010)
LCOR010 <- LCOR010[c(-1,-2,-8),]
LCOR010_fit = fitaci(LCOR010,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR010_fit)
summary(LCOR010_fit)
coef010 <- coef(LCOR010_fit)
tree_list <- append(tree_list, 'LCOR-010')
Vcmax_list <- append(Vcmax_list, coef010[1])
Jmax_list <- append(Jmax_list, coef010[2])
Rd_list <- append(Rd_list, coef010[3])


LCOR557 = read.xlsx(LC_spreadz,'LCOR-557')
viewdat(LCOR557)
LCOR557 <- LCOR557[c(-1,-2,-8),]
LCOR557_fit = fitaci(LCOR557,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR557_fit)
summary(LCOR557_fit)
coef557 <- coef(LCOR557_fit)
tree_list <- append(tree_list, 'LCOR-557')
Vcmax_list <- append(Vcmax_list, coef557[1])
Jmax_list <- append(Jmax_list, coef557[2])
Rd_list <- append(Rd_list, coef557[3])


LCOR197 = read.xlsx(LC_spreadz,'LCOR-197')
viewdat(LCOR197)
LCOR197 <- LCOR197[c(-1,-2,-8),]
LCOR197_fit = fitaci(LCOR197,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR197_fit)
summary(LCOR197_fit)
coef197 <- coef(LCOR197_fit)
tree_list <- append(tree_list, 'LCOR-197')
Vcmax_list <- append(Vcmax_list, coef197[1])
Jmax_list <- append(Jmax_list, coef197[2])
Rd_list <- append(Rd_list, coef197[3])


LCOR232 = read.xlsx(LC_spreadz,'LCOR-232')
viewdat(LCOR232)
LCOR232 <- LCOR232[c(-1,-2,-8),]
LCOR232_fit = fitaci(LCOR232,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
viewdat(LCOR232)
plotfits(LCOR232_fit)
summary(LCOR232_fit)
coef232 <- coef(LCOR232_fit)
coef232
tree_list <- append(tree_list, 'LCOR-232')
Vcmax_list <- append(Vcmax_list, coef232[1])
Jmax_list <- append(Jmax_list, coef232[2])
Rd_list <- append(Rd_list, coef232[3])


LCOR086 = read.xlsx(LC_spreadz,'LCOR-086')
viewdat(LCOR086)
LCOR086 <- LCOR086[c(-1,-2,-8),]
LCOR086_fit = fitaci(LCOR086,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR086_fit)
summary(LCOR086_fit)
coef086 <- coef(LCOR086_fit)
tree_list <- append(tree_list, 'LCOR-086')
Vcmax_list <- append(Vcmax_list, coef086[1])
Jmax_list <- append(Jmax_list, coef086[2])
Rd_list <- append(Rd_list, coef086[3])



LCOR511 = read.xlsx(LC_spreadz,'LCOR-511')
viewdat(LCOR511)
LCOR511 <- LCOR511[c(-1,-2,-8),]
LCOR511_fit = fitaci(LCOR511,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR511_fit)
summary(LCOR511_fit)
coef511 <- coef(LCOR511_fit)
tree_list <- append(tree_list, 'LCOR-511')
Vcmax_list <- append(Vcmax_list, coef511[1])
Jmax_list <- append(Jmax_list, coef511[2])
Rd_list <- append(Rd_list, coef511[3])


LCOR075 = read.xlsx(LC_spreadz,'LCOR-075')
viewdat(LCOR075)
LCOR075 <- LCOR075[c(-1,-2,-8),]
LCOR075_fit = fitaci(LCOR075,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR075_fit)
summary(LCOR075_fit)
coef075 <- coef(LCOR075_fit)
tree_list <- append(tree_list, 'LCOR-075')
Vcmax_list <- append(Vcmax_list, coef075[1])
Jmax_list <- append(Jmax_list, coef075[2])
Rd_list <- append(Rd_list, coef075[3])


LCOR081 = read.xlsx(LC_spreadz,'LCOR-081')
viewdat(LCOR081)
LCOR081 <- LCOR081[c(-1,-2,-8),]
LCOR081_fit = fitaci(LCOR081,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR081_fit)
summary(LCOR081_fit)
coef081 <- coef(LCOR081_fit)
tree_list <- append(tree_list, 'LCOR-081')
Vcmax_list <- append(Vcmax_list, coef081[1])
Jmax_list <- append(Jmax_list, coef081[2])
Rd_list <- append(Rd_list, coef081[3])


LCOR097 = read.xlsx(LC_spreadz,'LCOR-097')
viewdat(LCOR097)
LCOR097 <- LCOR097[c(-1,-2,-3,-8),]
LCOR097_fit = fitaci(LCOR097,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR097_fit)
summary(LCOR097_fit)
coef097 <- coef(LCOR097_fit)
#removed 300 value
tree_list <- append(tree_list, 'LCOR-097')
Vcmax_list <- append(Vcmax_list, coef097[1])
Jmax_list <- append(Jmax_list, coef097[2])
Rd_list <- append(Rd_list, coef097[3])

LCOR561 = read.xlsx(LC_spreadz,'LCOR-561')
viewdat(LCOR561)
LCOR561 <- LCOR561[c(-1,-2,-8),]
LCOR561_fit = fitaci(LCOR561,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR561_fit)
summary(LCOR561_fit)
coef561 <- coef(LCOR561_fit)
tree_list <- append(tree_list, 'LCOR-561')
Vcmax_list <- append(Vcmax_list, coef561[1])
Jmax_list <- append(Jmax_list, coef561[2])
Rd_list <- append(Rd_list, coef561[3])

LCOR615 = read.xlsx(LC_spreadz,'LCOR-615')
viewdat(LCOR615)
LCOR615 <- LCOR615[c(-1,-2,-8),]
LCOR615_fit = fitaci(LCOR615,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR615_fit)
summary(LCOR615_fit)
coef615 <- coef(LCOR615_fit)
tree_list <- append(tree_list, 'LCOR-615')
Vcmax_list <- append(Vcmax_list, coef615[1])
Jmax_list <- append(Jmax_list, coef615[2])
Rd_list <- append(Rd_list, coef615[3])

LCOR501 = read.xlsx(LC_spreadz,'LCOR-501')
viewdat(LCOR501)
LCOR501 <- LCOR501[c(-1,-2,-8),]
LCOR501_fit = fitaci(LCOR501,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR501_fit)
summary(LCOR501_fit)
coef501 <- coef(LCOR501_fit)
tree_list <- append(tree_list, 'LCOR-501')
Vcmax_list <- append(Vcmax_list, coef501[1])
Jmax_list <- append(Jmax_list, coef501[2])
Rd_list <- append(Rd_list, coef501[3])

LCOR286 = read.xlsx(LC_spreadz,'LCOR-286')
viewdat(LCOR286)
LCOR286 <- LCOR286[c(-1,-2,-8),]
LCOR286_fit = fitaci(LCOR286,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR286_fit)
summary(LCOR286_fit)
coef286 <- coef(LCOR286_fit)
tree_list <- append(tree_list, 'LCOR-286')
Vcmax_list <- append(Vcmax_list, coef286[1])
Jmax_list <- append(Jmax_list, coef286[2])
Rd_list <- append(Rd_list, coef286[3])

LCOR614 = read.xlsx(LC_spreadz,'LCOR-614')
viewdat(LCOR614)
LCOR614 <- LCOR614[c(-1,-2,-8),]
LCOR614_fit = fitaci(LCOR614,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR614_fit)
summary(LCOR614_fit)
coef614 <- coef(LCOR614_fit)
tree_list <- append(tree_list, 'LCOR-614')
Vcmax_list <- append(Vcmax_list, coef614[1])
Jmax_list <- append(Jmax_list, coef614[2])
Rd_list <- append(Rd_list, coef614[3])

LCOR495 = read.xlsx(LC_spreadz,'LCOR-495')
viewdat(LCOR495)
LCOR495 <- LCOR495[c(-1,-2,-8),]
LCOR495_fit = fitaci(LCOR495,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR495_fit)
summary(LCOR495_fit)
coef495 <- coef(LCOR495_fit)
tree_list <- append(tree_list, 'LCOR-495')
Vcmax_list <- append(Vcmax_list, coef495[1])
Jmax_list <- append(Jmax_list, coef495[2])
Rd_list <- append(Rd_list, coef495[3])

LCOR156 = read.xlsx(LC_spreadz,'LCOR-156')
viewdat(LCOR156)
LCOR156 <- LCOR156[c(-1,-2,-8),]
LCOR156_fit = fitaci(LCOR156,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR156_fit)
summary(LCOR156_fit)
coef156 <- coef(LCOR156_fit)
tree_list <- append(tree_list, 'LCOR-156')
Vcmax_list <- append(Vcmax_list, coef156[1])
Jmax_list <- append(Jmax_list, coef156[2])
Rd_list <- append(Rd_list, coef156[3])

LCOR031 = read.xlsx(LC_spreadz,'LCOR-031')
viewdat(LCOR031)
LCOR031 <- LCOR031[c(-1,-2,-8),]
LCOR031_fit = fitaci(LCOR031,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR031_fit)
summary(LCOR031_fit)
coef031 <- coef(LCOR031_fit)
tree_list <- append(tree_list, 'LCOR-031')
Vcmax_list <- append(Vcmax_list, coef031[1])
Jmax_list <- append(Jmax_list, coef031[2])
Rd_list <- append(Rd_list, coef031[3])

LCOR390 = read.xlsx(LC_spreadz,'LCOR-390')
viewdat(LCOR390)
LCOR390 <- LCOR390[c(-1,-2,-8),]
LCOR390_fit = fitaci(LCOR390,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR390_fit)
summary(LCOR390_fit)
coef390 <- coef(LCOR390_fit)
tree_list <- append(tree_list, 'LCOR-390')
Vcmax_list <- append(Vcmax_list, coef390[1])
Jmax_list <- append(Jmax_list, coef390[2])
Rd_list <- append(Rd_list, coef390[3])

LCOR202 = read.xlsx(LC_spreadz,'LCOR-202')
viewdat(LCOR202)
LCOR202 <- LCOR202[c(-1,-2,-8),]
LCOR202_fit = fitaci(LCOR202,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR202_fit)
summary(LCOR202_fit)
coef202 <- coef(LCOR202_fit)
tree_list <- append(tree_list, 'LCOR-202')
Vcmax_list <- append(Vcmax_list, coef202[1])
Jmax_list <- append(Jmax_list, coef202[2])
Rd_list <- append(Rd_list, coef202[3])

LCOR029 = read.xlsx(LC_spreadz,'LCOR-029')
viewdat(LCOR029)
LCOR029 <- LCOR029[c(-1,-2,-8),]
LCOR029_fit = fitaci(LCOR029,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR029_fit)
summary(LCOR029_fit)
coef029 <- coef(LCOR029_fit)
tree_list <- append(tree_list, 'LCOR-029')
Vcmax_list <- append(Vcmax_list, coef029[1])
Jmax_list <- append(Jmax_list, coef029[2])
Rd_list <- append(Rd_list, coef029[3])

LCOR342 = read.xlsx(LC_spreadz,'LCOR-342')
viewdat(LCOR342)
LCOR342 <- LCOR342[c(-1,-2,-8),]
LCOR342_fit = fitaci(LCOR342,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR342_fit)
summary(LCOR342_fit)
coef342 <- coef(LCOR342_fit)
tree_list <- append(tree_list, 'LCOR-342')
Vcmax_list <- append(Vcmax_list, coef342[1])
Jmax_list <- append(Jmax_list, coef342[2])
Rd_list <- append(Rd_list, coef342[3])

LCOR419 = read.xlsx(LC_spreadz,'LCOR-419')
viewdat(LCOR419)
LCOR419 <- LCOR419[c(-1,-2,-8),]
LCOR419_fit = fitaci(LCOR419,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR419_fit)
summary(LCOR419_fit)
coef419 <- coef(LCOR419_fit)
tree_list <- append(tree_list, 'LCOR-419')
Vcmax_list <- append(Vcmax_list, coef419[1])
Jmax_list <- append(Jmax_list, coef419[2])
Rd_list <- append(Rd_list, coef419[3])

LCOR463 = read.xlsx(LC_spreadz,'LCOR-463')
viewdat(LCOR463)
LCOR463 <- LCOR463[c(-1,-2,-8),]
LCOR463_fit = fitaci(LCOR463,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR463_fit)
summary(LCOR463_fit)
coef463 <- coef(LCOR463_fit)
tree_list <- append(tree_list, 'LCOR-463')
Vcmax_list <- append(Vcmax_list, coef463[1])
Jmax_list <- append(Jmax_list, coef463[2])
Rd_list <- append(Rd_list, coef463[3])

LCOR460 = read.xlsx(LC_spreadz,'LCOR-460')
viewdat(LCOR460)
LCOR460 <- LCOR460[c(-1,-2,-8),]
LCOR460_fit = fitaci(LCOR460,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR460_fit)
summary(LCOR460_fit)
coef460 <- coef(LCOR460_fit)
tree_list <- append(tree_list, 'LCOR-460')
Vcmax_list <- append(Vcmax_list, coef460[1])
Jmax_list <- append(Jmax_list, coef460[2])
Rd_list <- append(Rd_list, coef460[3])

LCOR577 = read.xlsx(LC_spreadz,'LCOR-577')
viewdat(LCOR577)
LCOR577 <- LCOR577[c(-1,-2,-8),]
LCOR577_fit = fitaci(LCOR577,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR577_fit)
summary(LCOR577_fit)
coef577 <- coef(LCOR577_fit)
tree_list <- append(tree_list, 'LCOR-577')
Vcmax_list <- append(Vcmax_list, coef577[1])
Jmax_list <- append(Jmax_list, coef577[2])
Rd_list <- append(Rd_list, coef577[3])

LCOR012 = read.xlsx(LC_spreadz,'LCOR-012')
viewdat(LCOR012)
LCOR012 <- LCOR012[c(-1,-2,-8),]
LCOR012_fit = fitaci(LCOR012,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR012_fit)
summary(LCOR012_fit)
coef012 <- coef(LCOR012_fit)
tree_list <- append(tree_list, 'LCOR-012')
Vcmax_list <- append(Vcmax_list, coef012[1])
Jmax_list <- append(Jmax_list, coef012[2])
Rd_list <- append(Rd_list, coef012[3])

LCOR611 = read.xlsx(LC_spreadz,'LCOR-611')
viewdat(LCOR611)
LCOR611 <- LCOR611[c(-1,-2,-8),]
LCOR611_fit = fitaci(LCOR611,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR611_fit)
summary(LCOR611_fit)
coef611 <- coef(LCOR611_fit)
tree_list <- append(tree_list, 'LCOR-611')
Vcmax_list <- append(Vcmax_list, coef611[1])
Jmax_list <- append(Jmax_list, coef611[2])
Rd_list <- append(Rd_list, coef611[3])

LCOR573 = read.xlsx(LC_spreadz,'LCOR-573')
viewdat(LCOR573)
LCOR573 <- LCOR573[c(-1,-2,-8),]
LCOR573_fit = fitaci(LCOR573,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR573_fit)
summary(LCOR573_fit)
coef573 <- coef(LCOR573_fit)
tree_list <- append(tree_list, 'LCOR-573')
Vcmax_list <- append(Vcmax_list, coef573[1])
Jmax_list <- append(Jmax_list, coef573[2])
Rd_list <- append(Rd_list, coef573[3])

LCOR148 = read.xlsx(LC_spreadz,'LCOR-148')
viewdat(LCOR148)
LCOR148 <- LCOR148[c(-1,-2,-8),]
LCOR148_fit = fitaci(LCOR148,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR148_fit)
summary(LCOR148_fit)
coef148 <- coef(LCOR148_fit)
tree_list <- append(tree_list, 'LCOR-148')
Vcmax_list <- append(Vcmax_list, coef148[1])
Jmax_list <- append(Jmax_list, coef148[2])
Rd_list <- append(Rd_list, coef148[3])

LCOR450 = read.xlsx(LC_spreadz,'LCOR-450')
viewdat(LCOR450)
LCOR450 <- LCOR450[c(-1,-8,-9,-14),]
LCOR450_fit = fitaci(LCOR450,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR450_fit)
summary(LCOR450_fit)
coef450 <- coef(LCOR450_fit)
tree_list <- append(tree_list, 'LCOR-450')
Vcmax_list <- append(Vcmax_list, coef450[1])
Jmax_list <- append(Jmax_list, coef450[2])
Rd_list <- append(Rd_list, coef450[3])

LCOR469 = read.xlsx(LC_spreadz,'LCOR-469')
viewdat(LCOR469)
LCOR469 <- LCOR469[c(-1,-2,-8,-14),]
LCOR469_fit = fitaci(LCOR469,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR469_fit)
summary(LCOR469_fit)
coef469 <- coef(LCOR469_fit)
tree_list <- append(tree_list, 'LCOR-469')
Vcmax_list <- append(Vcmax_list, coef469[1])
Jmax_list <- append(Jmax_list, coef469[2])
Rd_list <- append(Rd_list, coef469[3])

LCOR516 = read.xlsx(LC_spreadz,'LCOR-516')
viewdat(LCOR516)
LCOR516 <- LCOR516[c(-1,-2,-8),]
LCOR516_fit = fitaci(LCOR516,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR516_fit)
summary(LCOR516_fit)
coef516 <- coef(LCOR516_fit)
tree_list <- append(tree_list, 'LCOR-516')
Vcmax_list <- append(Vcmax_list, coef516[1])
Jmax_list <- append(Jmax_list, coef516[2])
Rd_list <- append(Rd_list, coef516[3])

LCOR411 = read.xlsx(LC_spreadz,'LCOR-411')
viewdat(LCOR411)
LCOR411 <- LCOR411[c(-1,-2,-8,-16),]
LCOR411_fit = fitaci(LCOR411,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR411_fit)
summary(LCOR411_fit)
coef411 <- coef(LCOR411_fit)
tree_list <- append(tree_list, 'LCOR-411')
Vcmax_list <- append(Vcmax_list, coef411[1])
Jmax_list <- append(Jmax_list, coef411[2])
Rd_list <- append(Rd_list, coef411[3])

LCOR209 = read.xlsx(LC_spreadz,'LCOR-209')
viewdat(LCOR209)
LCOR209 <- LCOR209[c(-1,-2,-8),]
LCOR209_fit = fitaci(LCOR209,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR209_fit)
summary(LCOR209_fit)
coef209 <- coef(LCOR209_fit)
tree_list <- append(tree_list, 'LCOR-209')
Vcmax_list <- append(Vcmax_list, coef209[1])
Jmax_list <- append(Jmax_list, coef209[2])
Rd_list <- append(Rd_list, coef209[3])

LCOR092 = read.xlsx(LC_spreadz,'LCOR-092')
viewdat(LCOR092)
LCOR092 <- LCOR092[c(-1,-2,-8,-14),]
LCOR092_fit = fitaci(LCOR092,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR092_fit)
summary(LCOR092_fit)
coef092 <- coef(LCOR092_fit)
tree_list <- append(tree_list, 'LCOR-092')
Vcmax_list <- append(Vcmax_list, coef092[1])
Jmax_list <- append(Jmax_list, coef092[2])
Rd_list <- append(Rd_list, coef092[3])

LCOR322 = read.xlsx(LC_spreadz,'LCOR-322')
viewdat(LCOR322)
LCOR322 <- LCOR322[c(-1,-2,-8,-15),]
LCOR322_fit = fitaci(LCOR322,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin', Rd = 'Rd'))
plotfits(LCOR322_fit)
summary(LCOR322_fit)
coef322 <- coef(LCOR322_fit)
tree_list <- append(tree_list, 'LCOR-322')
Vcmax_list <- append(Vcmax_list, coef322[1])
Jmax_list <- append(Jmax_list, coef322[2])
Rd_list <- append(Rd_list, coef322[3])

#compile parameter lists and id list into one dataframe

coef_lists = list(tree_list,
             Vcmax_list,
             Jmax_list,
             Rd_list)

# convert list of lists into dataframe
# by column
print(as.data.frame(do.call(cbind, coef_lists)))

efficiency_pars <- as.data.frame(do.call(cbind, coef_lists))
str(efficiency_pars)
colnames(efficiency_pars) = c('tree', 'Vcmax','Jmax','Rd')

write.csv(file = efficiency_pars)

efficiency_pars_II <- efficiency_pars
efficiency_pars_II$tree <- unlist(efficiency_pars_II$tree)
efficiency_pars_II$Vcmax <- unlist(efficiency_pars_II$Vcmax)
efficiency_pars_II$Jmax <- unlist(efficiency_pars_II$Jmax)
efficiency_pars_II$Rd <- unlist(efficiency_pars_II$Rd)
#combine with other parameter list via inner join

efficiency_pars_II$Jmax_Vcmax_ratio <- efficiency_pars_II$Jmax/efficiency_pars_II$Vcmax

aci_parameters_list <- inner_join(efficiency_pars_II, all_trees_wide, by = "tree")
aci_parameters_list_flr <- inner_join(efficiency_pars_II, all_trees_flr_wide, by = "tree")


write.csv(aci_parameters_list, file = "2022_physiology_analysis/CO2_response/CO2_response_modeling/Aci_parameters_list.csv")
write.csv(aci_parameters_list_flr, file = "2022_physiology_analysis/CO2_response/CO2_response_modeling/Aci_parameters_list_flr.csv")

