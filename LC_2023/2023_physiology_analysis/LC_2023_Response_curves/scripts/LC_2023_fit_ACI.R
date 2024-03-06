#process fits for A-ci curves and output them in a csv

#install.packages("plantecophys")
library(plantecophys)
#install.packages("nlstools")
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
#install.packages("openxlsx")
library(openxlsx)
library(mgcv)
#install.packages("photosynthesis")
library(photosynthesis)
#read in compiled datasets

ACI_compiled <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/CO2_response_curves_compiled.csv")

#long format

#all_trees_wide <- pivot_wider(ACI_compiled,id_cols = c("ID"), names_from = Csetpoint, values_from = c(A,gsw,Ci))

#fx for viewing data
viewdat = function(dat){plot(y = dat$A, x = dat$Ci);abline(h = seq(0,20,1));dat}

#fx for plotting fits

plotfits = function(fit){
  jpeg(filename = paste('',deparse(substitute(fit)),'.jpeg',sep=''),width = 10,height=7,units='in',res = 100)
  plot(fit,main =deparse(substitute(fit)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel),labels = paste('Vcmax =',round(fit$pars[1,1],2)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel)-1,labels = paste('Jmax =',round(fit$pars[2,1],2)))
  dev.off()
}

#options("device")
#create lists to store data in

tree_list = list()
Vcmax_list = list()
Jmax_list = list()
Rd_list = list()
Amax_list = list()
Sample_date = list()

##look at fits, model and graph them and get outputs and assign to storage dataframe

LCOR307 = subset(ACI_compiled, ID == "LCOR-307")
LCOR307_1 = LCOR307[1:14,]
LCOR307_2 = LCOR307[15:28,]
viewdat(LCOR307_1)
viewdat(LCOR307_2)
#LCOR307 <- LCOR307[c(-1,-2,-7),]
LCOR307_fit_1 = fitaci(LCOR307_1,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
LCOR307_fit_2 = fitaci(LCOR307_2,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))

#fit_aci_response(LCOR307_1,varnames = list(A_net = "A",T_leaf = "Tleaf",C_i = "Ci",PPFD = "Qin",g_mc))
#?fit_aci_response
plotfits(LCOR307_fit_1)
plotfits(LCOR307_fit_2)
summary(LCOR307_fit_1)
plotfits(LCOR307_fit_2)
coef307_1 <- coef(LCOR307_fit_1)
coef307_2 <- coef(LCOR307_fit_2)
coef307_1
coef307_2
tree_list <- append(tree_list, 'LCOR-307')
Vcmax_list <- append(Vcmax_list, coef307_1[1])
Jmax_list <- append(Jmax_list, coef307_1[2])
Rd_list <- append(Rd_list, coef307_1[3])
Amax_list <- append(Amax_list, max(LCOR307_fit_1$df$Amodel))
tree_list <- append(tree_list, 'LCOR-307_2')
Vcmax_list <- append(Vcmax_list, coef307_2[1])
Jmax_list <- append(Jmax_list, coef307_2[2])
Rd_list <- append(Rd_list, coef307_2[3])
Amax_list <- append(Amax_list, max(LCOR307_fit_2$df$Amodel))

LCOR318 = subset(ACI_compiled, ID == "LCOR-318")
viewdat(LCOR318)
#LCOR318 <- LCOR318[c(-1,-2,-7),]
LCOR318_fit = fitaci(LCOR318,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR318_fit)
summary(LCOR318_fit)
coef318 <- coef(LCOR318_fit)
LCOR318_fit$Photosyn(Ci=200)
coef318
tree_list <- append(tree_list, 'LCOR-318')
Vcmax_list <- append(Vcmax_list, coef318[1])
Jmax_list <- append(Jmax_list, coef318[2])
Rd_list <- append(Rd_list, coef318[3])
Amax_list <- append(Amax_list, max(LCOR318_fit$df$Amodel))

LCOR412 = subset(ACI_compiled, ID == "LCOR-412")
viewdat(LCOR412)
#LCOR412 <- LCOR412[c(-1,-2,-7),]
LCOR412_fit = fitaci(LCOR412,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR412_fit)
summary(LCOR412_fit)
coef412 <- coef(LCOR412_fit)
LCOR412_fit$Photosyn(Ci=200)
coef412
tree_list <- append(tree_list, 'LCOR-412')
Vcmax_list <- append(Vcmax_list, coef412[1])
Jmax_list <- append(Jmax_list, coef412[2])
Rd_list <- append(Rd_list, coef412[3])
Amax_list <- append(Amax_list, max(LCOR412_fit$df$Amodel))

LCOR063 = subset(ACI_compiled, ID == "LCOR-063")
viewdat(LCOR063)
LCOR063 <- LCOR063[c(-11,-12),]
LCOR063_fit = fitaci(LCOR063,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR063_fit)
summary(LCOR063_fit)
coef063 <- coef(LCOR063_fit)
LCOR063_fit$Photosyn(Ci=200)
coef063
tree_list <- append(tree_list, 'LCOR-063')
Vcmax_list <- append(Vcmax_list, coef063[1])
Jmax_list <- append(Jmax_list, coef063[2])
Rd_list <- append(Rd_list, coef063[3])
Amax_list <- append(Amax_list, max(LCOR063_fit$df$Amodel))

LCOR546 = subset(ACI_compiled, ID == "LCOR-546")
viewdat(LCOR546)
LCOR546_fit = fitaci(LCOR546,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR546_fit)
summary(LCOR546_fit)
coef546 <- coef(LCOR546_fit)
LCOR546_fit$Photosyn(Ci=200)
coef546
tree_list <- append(tree_list, 'LCOR-546')
Vcmax_list <- append(Vcmax_list, coef546[1])
Jmax_list <- append(Jmax_list, coef546[2])
Rd_list <- append(Rd_list, coef546[3])
Amax_list <- append(Amax_list, max(LCOR546_fit$df$Amodel))

LCOR347 = subset(ACI_compiled, ID == "LCOR-347")
viewdat(LCOR347)
LCOR347_fit = fitaci(LCOR347,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR347_fit)
summary(LCOR347_fit)
coef347 <- coef(LCOR347_fit)
LCOR347_fit$Photosyn(Ci=200)
coef347
tree_list <- append(tree_list, 'LCOR-347')
Vcmax_list <- append(Vcmax_list, coef347[1])
Jmax_list <- append(Jmax_list, coef347[2])
Rd_list <- append(Rd_list, coef347[3])
Amax_list <- append(Amax_list, max(LCOR347_fit$df$Amodel))

LCOR268 = subset(ACI_compiled, ID == "LCOR-268")
viewdat(LCOR268)
LCOR268_fit = fitaci(LCOR268,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR268_fit)
summary(LCOR268_fit)
coef268 <- coef(LCOR268_fit)
LCOR268_fit$Photosyn(Ci=200)
coef268
tree_list <- append(tree_list, 'LCOR-268')
Vcmax_list <- append(Vcmax_list, coef268[1])
Jmax_list <- append(Jmax_list, coef268[2])
Rd_list <- append(Rd_list, coef268[3])
Amax_list <- append(Amax_list, max(LCOR268_fit$df$Amodel))


LCOR155 = subset(ACI_compiled, ID == "LCOR-155")
viewdat(LCOR155)

LCOR155_fit = fitaci(LCOR155,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR155_fit)
summary(LCOR155_fit)
coef155 <- coef(LCOR155_fit)
LCOR155_fit$Photosyn(Ci=200)
coef155
tree_list <- append(tree_list, 'LCOR-155')
Vcmax_list <- append(Vcmax_list, coef155[1])
Jmax_list <- append(Jmax_list, coef155[2])
Rd_list <- append(Rd_list, coef155[3])
Amax_list <- append(Amax_list, max(LCOR155_fit$df$Amodel))

LCOR322 = subset(ACI_compiled, ID == "LCOR-322")
viewdat(LCOR322)

LCOR322_fit = fitaci(LCOR322,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR322_fit)
summary(LCOR322_fit)
coef322 <- coef(LCOR322_fit)
LCOR322_fit$Photosyn(Ci=200)
coef322
tree_list <- append(tree_list, 'LCOR-322')
Vcmax_list <- append(Vcmax_list, coef322[1])
Jmax_list <- append(Jmax_list, coef322[2])
Rd_list <- append(Rd_list, coef322[3])
Amax_list <- append(Amax_list, max(LCOR322_fit$df$Amodel))

LCOR073 = subset(ACI_compiled, ID == "LCOR-073")
viewdat(LCOR073)

LCOR073_fit = fitaci(LCOR073,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR073_fit)
summary(LCOR073_fit)
coef073 <- coef(LCOR073_fit)
LCOR073_fit$Photosyn(Ci=200)
coef073
tree_list <- append(tree_list, 'LCOR-073')
Vcmax_list <- append(Vcmax_list, coef073[1])
Jmax_list <- append(Jmax_list, coef073[2])
Rd_list <- append(Rd_list, coef073[3])
Amax_list <- append(Amax_list, max(LCOR073_fit$df$Amodel))

LCOR164 = subset(ACI_compiled, ID == "LCOR-164")
viewdat(LCOR164)

LCOR164 = subset(ACI_compiled, ID == "LCOR-164")
LCOR164_1 = LCOR164[1:14,]
LCOR164_2 = LCOR164[15:28,]
viewdat(LCOR164_1)
viewdat(LCOR164_2)
#LCOR307 <- LCOR307[c(-1,-2,-7),]
LCOR164_fit_1 = fitaci(LCOR164_1,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
LCOR164_fit_2 = fitaci(LCOR164_2,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR164_fit_1)
plotfits(LCOR164_fit_2)
summary(LCOR164_fit_1)
summary(LCOR164_fit_2)
coef164_1 <- coef(LCOR164_fit_1)
coef164_2 <- coef(LCOR164_fit_2)
coef164_1
coef164_2
tree_list <- append(tree_list, 'LCOR-164_1')
Vcmax_list <- append(Vcmax_list, coef164_1[1])
Jmax_list <- append(Jmax_list, coef164_1[2])
Rd_list <- append(Rd_list, coef164_1[3])
Amax_list <- append(Amax_list, max(LCOR164_fit_1$df$Amodel))
tree_list <- append(tree_list, 'LCOR-164_2')
Vcmax_list <- append(Vcmax_list, coef164_2[1])
Jmax_list <- append(Jmax_list, coef164_2[2])
Rd_list <- append(Rd_list, coef164_2[3])
Amax_list <- append(Amax_list, max(LCOR164_fit_2$df$Amodel))

LCOR199 = subset(ACI_compiled, ID == "LCOR-199")
viewdat(LCOR199)

LCOR199_fit = fitaci(LCOR199,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR199_fit)
summary(LCOR199_fit)
coef199 <- coef(LCOR199_fit)
LCOR199_fit$Photosyn(Ci=200)
coef199
tree_list <- append(tree_list, 'LCOR-199')
Vcmax_list <- append(Vcmax_list, coef199[1])
Jmax_list <- append(Jmax_list, coef199[2])
Rd_list <- append(Rd_list, coef199[3])
Amax_list <- append(Amax_list, max(LCOR199_fit$df$Amodel))

LCOR573 = subset(ACI_compiled, ID == "LCOR-573")
viewdat(LCOR573)

LCOR573_fit = fitaci(LCOR573,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR573_fit)
summary(LCOR573_fit)
coef573 <- coef(LCOR573_fit)
LCOR573_fit$Photosyn(Ci=200)
coef573
tree_list <- append(tree_list, 'LCOR-573')
Vcmax_list <- append(Vcmax_list, coef573[1])
Jmax_list <- append(Jmax_list, coef573[2])
Rd_list <- append(Rd_list, coef573[3])
Amax_list <- append(Amax_list, max(LCOR573_fit$df$Amodel))


LCOR456 = subset(ACI_compiled, ID == "LCOR-456")
viewdat(LCOR456)

LCOR456_fit = fitaci(LCOR456,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR456_fit)
summary(LCOR456_fit)
coef456 <- coef(LCOR456_fit)
LCOR456_fit$Photosyn(Ci=200)
coef456
tree_list <- append(tree_list, 'LCOR-456')
Vcmax_list <- append(Vcmax_list, coef456[1])
Jmax_list <- append(Jmax_list, coef456[2])
Rd_list <- append(Rd_list, coef456[3])
Amax_list <- append(Amax_list, max(LCOR456_fit$df$Amodel))

LCOR089 = subset(ACI_compiled, ID == "LCOR-089")
viewdat(LCOR089)
#does not look good!

LCOR089_fit = fitaci(LCOR089,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR089_fit)
summary(LCOR089_fit)
coef089 <- coef(LCOR089_fit)
LCOR089_fit$Photosyn(Ci=200)
coef089
#tree_list <- append(tree_list, 'LCOR-089')
#Vcmax_list <- append(Vcmax_list, coef089[1])
#Jmax_list <- append(Jmax_list, coef089[2])
#Rd_list <- append(Rd_list, coef089[3])

LCOR229 = subset(ACI_compiled, ID == "LCOR-229")
viewdat(LCOR229)
#LCOR229 <- LCOR229[c(-11),]
LCOR229_fit = fitaci(LCOR229,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR229_fit)
summary(LCOR229_fit)
coef229 <- coef(LCOR229_fit)
LCOR229_fit$Photosyn(Ci=200)
coef229
tree_list <- append(tree_list, 'LCOR-229')
Vcmax_list <- append(Vcmax_list, coef229[1])
Jmax_list <- append(Jmax_list, coef229[2])
Rd_list <- append(Rd_list, coef229[3])
Amax_list <- append(Amax_list, max(LCOR229_fit$df$Amodel))

LCOR228 = subset(ACI_compiled, ID == "LCOR-228")
LCOR228 <- LCOR228[c(-5,-1,-4,-7),]
viewdat(LCOR228)

LCOR228_fit = fitaci(LCOR228,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR228_fit)
summary(LCOR228_fit)
coef228 <- coef(LCOR228_fit)
LCOR228_fit$Photosyn(Ci=200)
coef228
tree_list <- append(tree_list, 'LCOR-228')
Vcmax_list <- append(Vcmax_list, coef228[1])
Jmax_list <- append(Jmax_list, coef228[2])
Rd_list <- append(Rd_list, coef228[3])
Amax_list <- append(Amax_list, max(LCOR228_fit$df$Amodel))

LCOR577 = subset(ACI_compiled, ID == "LCOR-577")
LCOR577 <- LCOR577[c(-8,-7,-6,-1),]
viewdat(LCOR577)
LCOR577_fit = fitaci(LCOR577,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR577_fit)
summary(LCOR577_fit)
coef577 <- coef(LCOR577_fit)
LCOR577_fit$Photosyn(Ci=200)
coef577
tree_list <- append(tree_list, 'LCOR-577')
Vcmax_list <- append(Vcmax_list, coef577[1])
Jmax_list <- append(Jmax_list, coef577[2])
Rd_list <- append(Rd_list, coef577[3])
Amax_list <- append(Amax_list, max(LCOR577_fit$df$Amodel))

LCOR611 = subset(ACI_compiled, ID == "LCOR-611")
LCOR611 <- LCOR611[c(-1,-12),]
viewdat(LCOR611)
LCOR611_fit = fitaci(LCOR611,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR611_fit)
summary(LCOR611_fit)
coef611 <- coef(LCOR611_fit)
LCOR611_fit$Photosyn(Ci=200)
coef611
tree_list <- append(tree_list, 'LCOR-611')
Vcmax_list <- append(Vcmax_list, coef611[1])
Jmax_list <- append(Jmax_list, coef611[2])
Rd_list <- append(Rd_list, coef611[3])
Amax_list <- append(Amax_list, max(LCOR611_fit$df$Amodel))

LCOR280 = subset(ACI_compiled, ID == "LCOR-280")

viewdat(LCOR280)
LCOR280_fit = fitaci(LCOR280,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR280_fit)
summary(LCOR280_fit)
coef280 <- coef(LCOR280_fit)
LCOR280_fit$Photosyn(Ci=200)
coef280
tree_list <- append(tree_list, 'LCOR-280')
Vcmax_list <- append(Vcmax_list, coef280[1])
Jmax_list <- append(Jmax_list, coef280[2])
Rd_list <- append(Rd_list, coef280[3])
Amax_list <- append(Amax_list, max(LCOR280_fit$df$Amodel))


LCOR291 = subset(ACI_compiled, ID == "LCOR-291")
LCOR291_1 = LCOR291[1:14,]
LCOR291_1 <- LCOR291_1[c(-1),]
viewdat(LCOR291_1)
LCOR291_2 = LCOR291[15:28,]
viewdat(LCOR291_2)
#LCOR307 <- LCOR307[c(-1,-2,-7),]
LCOR291_fit_1 = fitaci(LCOR291_1,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
LCOR291_fit_2 = fitaci(LCOR291_2,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR291_fit_1)
plotfits(LCOR291_fit_2)
summary(LCOR291_fit_1)
summary(LCOR291_fit_2)
coef291_1 <- coef(LCOR291_fit_1)
coef291_2 <- coef(LCOR291_fit_2)
coef291_1
coef291_2
tree_list <- append(tree_list, 'LCOR-291')
Vcmax_list <- append(Vcmax_list, coef291_1[1])
Jmax_list <- append(Jmax_list, coef291_1[2])
Rd_list <- append(Rd_list, coef291_1[3])
tree_list <- append(tree_list, 'LCOR-291_2')
Vcmax_list <- append(Vcmax_list, coef291_2[1])
Jmax_list <- append(Jmax_list, coef291_2[2])
Rd_list <- append(Rd_list, coef291_2[3])
Amax_list <- append(Amax_list, max(LCOR291_fit_1$df$Amodel))
Amax_list <- append(Amax_list, max(LCOR291_fit_2$df$Amodel))

LCOR506 = subset(ACI_compiled, ID == "LCOR-506")
LCOR506 <- LCOR506[c(-1,-2),]
viewdat(LCOR506)
#looks bad
LCOR506_fit = fitaci(LCOR506,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR506_fit)
summary(LCOR506_fit)
coef506 <- coef(LCOR506_fit)
LCOR506_fit$Photosyn(Ci=200)
coef506
#tree_list <- append(tree_list, 'LCOR-506')
#Vcmax_list <- append(Vcmax_list, coef506[1])
#Jmax_list <- append(Jmax_list, coef506[2])
#Rd_list <- append(Rd_list, coef506[3])

LCOR501 = subset(ACI_compiled, ID == "LCOR-501")
LCOR501 <- LCOR501[c(-1,-2),]
viewdat(LCOR501)
LCOR501_fit = fitaci(LCOR501,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR501_fit)
summary(LCOR501_fit)
coef501 <- coef(LCOR501_fit)
LCOR501_fit$Photosyn(Ci=200)
coef501
tree_list <- append(tree_list, 'LCOR-501')
Vcmax_list <- append(Vcmax_list, coef501[1])
Jmax_list <- append(Jmax_list, coef501[2])
Rd_list <- append(Rd_list, coef501[3])
Amax_list <- append(Amax_list, max(LCOR501_fit$df$Amodel))

LCOR211 = subset(ACI_compiled, ID == "LCOR-211")
LCOR211 <- LCOR211[c(-7,-9),]
viewdat(LCOR211)
LCOR211_fit = fitaci(LCOR211,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR211_fit)
summary(LCOR211_fit)
coef211 <- coef(LCOR211_fit)
LCOR211_fit$Photosyn(Ci=200)
coef211
tree_list <- append(tree_list, 'LCOR-211')
Vcmax_list <- append(Vcmax_list, coef211[1])
Jmax_list <- append(Jmax_list, coef211[2])
Rd_list <- append(Rd_list, coef211[3])
Amax_list <- append(Amax_list, max(LCOR211_fit$df$Amodel))

LCOR274 = subset(ACI_compiled, ID == "LCOR-274")
LCOR274 <- LCOR274[c(-7),]
viewdat(LCOR274)
#looks C4 like
LCOR274_fit = fitaci(LCOR274,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR274_fit)
summary(LCOR274_fit)
coef274 <- coef(LCOR274_fit)
LCOR274_fit$Photosyn(Ci=200)
coef274
tree_list <- append(tree_list, 'LCOR-274')
Vcmax_list <- append(Vcmax_list, coef274[1])
Jmax_list <- append(Jmax_list, coef274[2])
Rd_list <- append(Rd_list, coef274[3])
Amax_list <- append(Amax_list, max(LCOR274_fit$df$Amodel))

LCOR231 = subset(ACI_compiled, ID == "LCOR-231")
LCOR231 <- LCOR231[c(-7),]
viewdat(LCOR231)
LCOR231_fit = fitaci(LCOR231,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR231_fit)
summary(LCOR231_fit)
coef231 <- coef(LCOR231_fit)
LCOR231_fit$Photosyn(Ci=200)
coef231
tree_list <- append(tree_list, 'LCOR-231')
Vcmax_list <- append(Vcmax_list, coef231[1])
Jmax_list <- append(Jmax_list, coef231[2])
Rd_list <- append(Rd_list, coef231[3])
Amax_list <- append(Amax_list, max(LCOR231_fit$df$Amodel))

LCOR101 = subset(ACI_compiled, ID == "LCOR-101")
LCOR101 <- LCOR101[c(-3),]
viewdat(LCOR101)
LCOR101_fit = fitaci(LCOR101,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR101_fit)
summary(LCOR101_fit)
coef101 <- coef(LCOR101_fit)
LCOR101_fit$Photosyn(Ci=200)
coef101
tree_list <- append(tree_list, 'LCOR-101')
Vcmax_list <- append(Vcmax_list, coef101[1])
Jmax_list <- append(Jmax_list, coef101[2])
Rd_list <- append(Rd_list, coef101[3])
Amax_list <- append(Amax_list, max(LCOR101_fit$df$Amodel))


LCOR204 = subset(ACI_compiled, ID == "LCOR-204")
LCOR204 <- LCOR204[c(-15),]
viewdat(LCOR204)
LCOR204_fit = fitaci(LCOR204,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR204_fit)
summary(LCOR204_fit)
coef204 <- coef(LCOR204_fit)
LCOR204_fit$Photosyn(Ci=200)
coef204
tree_list <- append(tree_list, 'LCOR-204')
Vcmax_list <- append(Vcmax_list, coef204[1])
Jmax_list <- append(Jmax_list, coef204[2])
Rd_list <- append(Rd_list, coef204[3])
Amax_list <- append(Amax_list, max(LCOR204_fit$df$Amodel))

LCOR203 = subset(ACI_compiled, ID == "LCOR-203")

viewdat(LCOR203)
LCOR203_fit = fitaci(LCOR203,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR203_fit)
summary(LCOR203_fit)
coef203 <- coef(LCOR203_fit)
LCOR203_fit$Photosyn(Ci=200)
coef203
tree_list <- append(tree_list, 'LCOR-203')
Vcmax_list <- append(Vcmax_list, coef203[1])
Jmax_list <- append(Jmax_list, coef203[2])
Rd_list <- append(Rd_list, coef203[3])
Amax_list <- append(Amax_list, max(LCOR203_fit$df$Amodel))

LCOR610 = subset(ACI_compiled, ID == "LCOR-610")
viewdat(LCOR610)
LCOR610_fit = fitaci(LCOR610,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR610_fit)
summary(LCOR610_fit)
coef610 <- coef(LCOR610_fit)
LCOR610_fit$Photosyn(Ci=200)
coef610
tree_list <- append(tree_list, 'LCOR-610')
Vcmax_list <- append(Vcmax_list, coef610[1])
Jmax_list <- append(Jmax_list, coef610[2])
Rd_list <- append(Rd_list, coef610[3])
Amax_list <- append(Amax_list, max(LCOR610_fit$df$Amodel))

LCOR275 = subset(ACI_compiled, ID == "LCOR-275")

viewdat(LCOR275)
LCOR275_fit = fitaci(LCOR275,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR275_fit)
summary(LCOR275_fit)
coef275 <- coef(LCOR275_fit)
LCOR275_fit$Photosyn(Ci=200)
coef275
tree_list <- append(tree_list, 'LCOR-275')
Vcmax_list <- append(Vcmax_list, coef275[1])
Jmax_list <- append(Jmax_list, coef275[2])
Rd_list <- append(Rd_list, coef275[3])
Amax_list <- append(Amax_list, max(LCOR275_fit$df$Amodel))

LCOR516 = subset(ACI_compiled, ID == "LCOR-516")
LCOR516 <- LCOR516[c(-11,-1),]
viewdat(LCOR516)
#does not look good
LCOR516_fit = fitaci(LCOR516,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR516_fit)
summary(LCOR516_fit)
coef516 <- coef(LCOR516_fit)
LCOR516_fit$Photosyn(Ci=200)
coef516
#tree_list <- append(tree_list, 'LCOR-516')
#Vcmax_list <- append(Vcmax_list, coef516[1])
#Jmax_list <- append(Jmax_list, coef516[2])
#Rd_list <- append(Rd_list, coef516[3])

LCOR003 = subset(ACI_compiled, ID == "LCOR-003")
LCOR003 <- LCOR003[c(-1,-7,-12),]
viewdat(LCOR003)
LCOR003_fit = fitaci(LCOR003,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR003_fit)
summary(LCOR003_fit)
coef003 <- coef(LCOR003_fit)
LCOR003_fit$Photosyn(Ci=200)
coef003
tree_list <- append(tree_list, 'LCOR-003')
Vcmax_list <- append(Vcmax_list, coef003[1])
Jmax_list <- append(Jmax_list, coef003[2])
Rd_list <- append(Rd_list, coef003[3])
Amax_list <- append(Amax_list, max(LCOR003_fit$df$Amodel))


LCOR578 = subset(ACI_compiled, ID == "LCOR-578")
LCOR578 <- LCOR578[c(-1,-12),]
viewdat(LCOR578)
LCOR578_fit = fitaci(LCOR578,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR578_fit)
summary(LCOR578_fit)
coef578 <- coef(LCOR578_fit)
LCOR578_fit$Photosyn(Ci=200)
coef578
tree_list <- append(tree_list, 'LCOR-578')
Vcmax_list <- append(Vcmax_list, coef578[1])
Jmax_list <- append(Jmax_list, coef578[2])
Rd_list <- append(Rd_list, coef578[3])
Amax_list <- append(Amax_list, max(LCOR578_fit$df$Amodel))

LCOR417 = subset(ACI_compiled, ID == "LCOR-417")
LCOR417 <- LCOR417[c(-8),]
viewdat(LCOR417)
#does not look good
LCOR417_fit = fitaci(LCOR417,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR417_fit)
summary(LCOR417_fit)
coef417 <- coef(LCOR417_fit)
LCOR417_fit$Photosyn(Ci=200)
coef417
#strange
#tree_list <- append(tree_list, 'LCOR-417')
#Vcmax_list <- append(Vcmax_list, coef417[1])
#Jmax_list <- append(Jmax_list, coef417[2])
#Rd_list <- append(Rd_list, coef417[3])

LCOR303 = subset(ACI_compiled, ID == "LCOR-303")
LCOR303 <- LCOR303[c(-1,-2),]
viewdat(LCOR303)
LCOR303_fit = fitaci(LCOR303,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR303_fit)
summary(LCOR303_fit)
coef303 <- coef(LCOR303_fit)
LCOR303_fit$Photosyn(Ci=200)
coef303
tree_list <- append(tree_list, 'LCOR-303')
Vcmax_list <- append(Vcmax_list, coef303[1])
Jmax_list <- append(Jmax_list, coef303[2])
Rd_list <- append(Rd_list, coef303[3])
Amax_list <- append(Amax_list, max(LCOR303_fit$df$Amodel))

LCOR502 = subset(ACI_compiled, ID == "LCOR-502")
LCOR502 <- LCOR502[c(-2),]
viewdat(LCOR502)
#does nooottt look good
#LCOR502_fit = fitaci(LCOR502,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
#plotfits(LCOR502_fit)
#summary(LCOR502_fit)
#coef502 <- coef(LCOR502_fit)
#LCOR502_fit$Photosyn(Ci=200)
#coef502
#tree_list <- append(tree_list, 'LCOR-502')
#Vcmax_list <- append(Vcmax_list, coef502[1])
#Jmax_list <- append(Jmax_list, coef502[2])
#Rd_list <- append(Rd_list, coef502[3])

LCOR314 = subset(ACI_compiled, ID == "LCOR-314")

viewdat(LCOR314)
LCOR314_fit = fitaci(LCOR314,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR314_fit)
summary(LCOR314_fit)
coef314 <- coef(LCOR314_fit)
LCOR314_fit$Photosyn(Ci=200)
coef314
tree_list <- append(tree_list, 'LCOR-314')
Vcmax_list <- append(Vcmax_list, coef314[1])
Jmax_list <- append(Jmax_list, coef314[2])
Rd_list <- append(Rd_list, coef314[3])
Amax_list <- append(Amax_list, max(LCOR314_fit$df$Amodel))

LCOR092 = subset(ACI_compiled, ID == "LCOR-092")

viewdat(LCOR092)
LCOR092_fit = fitaci(LCOR092,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR092_fit)
summary(LCOR092_fit)
coef092 <- coef(LCOR092_fit)
LCOR092_fit$Photosyn(Ci=200)
coef092
tree_list <- append(tree_list, 'LCOR-092')
Vcmax_list <- append(Vcmax_list, coef092[1])
Jmax_list <- append(Jmax_list, coef092[2])
Rd_list <- append(Rd_list, coef092[3])
Amax_list <- append(Amax_list, max(LCOR092_fit$df$Amodel))

LCOR065 = subset(ACI_compiled, ID == "LCOR-065")
LCOR065 <- LCOR065[c(-2,-1),]
viewdat(LCOR065)
LCOR065_fit = fitaci(LCOR065,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR065_fit)
summary(LCOR065_fit)
coef065 <- coef(LCOR065_fit)
LCOR065_fit$Photosyn(Ci=200)
coef065
tree_list <- append(tree_list, 'LCOR-065')
Vcmax_list <- append(Vcmax_list, coef065[1])
Jmax_list <- append(Jmax_list, coef065[2])
Rd_list <- append(Rd_list, coef065[3])
Amax_list <- append(Amax_list, max(LCOR065_fit$df$Amodel))

LCOR165 = subset(ACI_compiled, ID == "LCOR-165")
LCOR165 <- LCOR165[c(-8,-2),]
viewdat(LCOR165)
LCOR165_fit = fitaci(LCOR165,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR165_fit)
summary(LCOR165_fit)
coef165 <- coef(LCOR165_fit)
LCOR165_fit$Photosyn(Ci=200)
coef165
tree_list <- append(tree_list, 'LCOR-165')
Vcmax_list <- append(Vcmax_list, coef165[1])
Jmax_list <- append(Jmax_list, coef165[2])
Rd_list <- append(Rd_list, coef165[3])
Amax_list <- append(Amax_list, max(LCOR165_fit$df$Amodel))

LCOR242 = subset(ACI_compiled, ID == "LCOR-242")
viewdat(LCOR242)
LCOR242_fit = fitaci(LCOR242,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR242_fit)
summary(LCOR242_fit)
coef242 <- coef(LCOR242_fit)
LCOR242_fit$Photosyn(Ci=200)
coef242
tree_list <- append(tree_list, 'LCOR-242')
Vcmax_list <- append(Vcmax_list, coef242[1])
Jmax_list <- append(Jmax_list, coef242[2])
Rd_list <- append(Rd_list, coef242[3])
Amax_list <- append(Amax_list, max(LCOR242_fit$df$Amodel))

LCOR040 = subset(ACI_compiled, ID == "LCOR-040")

viewdat(LCOR040)
LCOR040_fit = fitaci(LCOR040,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
plotfits(LCOR040_fit)
summary(LCOR040_fit)
coef040 <- coef(LCOR040_fit)
LCOR040_fit$Photosyn(Ci=200)
coef040
tree_list <- append(tree_list, 'LCOR-040')
Vcmax_list <- append(Vcmax_list, coef040[1])
Jmax_list <- append(Jmax_list, coef040[2])
Rd_list <- append(Rd_list, coef040[3])
Amax_list <- append(Amax_list, max(LCOR040_fit$df$Amodel))




#compile parameter lists and id list into one dataframe
tree_list <- unlist(tree_list)
tree_list <- as.data.frame(as.matrix(tree_list))
Vcmax_list <- unlist(Vcmax_list)
Vcmax_list <- as.data.frame(as.matrix(Vcmax_list))
Jmax_list <- unlist(Jmax_list)
Jmax_list <- as.data.frame(as.matrix(Jmax_list))
Rd_list <- unlist(Rd_list)
Rd_list <- as.data.frame(as.matrix(Rd_list))
Amax_list <- unlist(Amax_list)
Amax_list <- as.data.frame(as.matrix(Amax_list))

efficiency_pars <- cbind(tree_list,Vcmax_list,Jmax_list,Rd_list,Amax_list)

colnames(efficiency_pars) = c('ID', 'Vcmax','Jmax','Rd',"Amax")
efficiency_pars$ID_2 <- substr(efficiency_pars$ID,1,8)
str(efficiency_pars)

efficiency_pars2 <- efficiency_pars %>% group_by(ID_2) %>% summarize(
  Vcmax = mean(Vcmax),
  Jmax = mean(Jmax),
  Rd = mean(Rd),
  Amax = mean(Amax)
)

efficiency_pars2$ID <- efficiency_pars2$ID_2

efficiency_pars2$Jmax_Vcmax_ratio <- efficiency_pars2$Jmax/efficiency_pars2$Vcmax

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2"))


aci_parameters_list <- inner_join(efficiency_pars2, LC_meta, by = "ID")
sample_date <- subset(ACI_compiled, select = c("ID","sample_date"))

aci_parameters_list2 <- distinct(inner_join(aci_parameters_list, sample_date, by = "ID"))

aci_parameters_event_summary <- aci_parameters_list2 %>% group_by(event_short) %>% summarize(
  n = n(),
  Vcmax_sd = sd(Vcmax),
  Jmax_sd = sd(Jmax),
  Vcmax = mean(Vcmax),
  Jmax = mean(Jmax),
  Vcmax_se = Vcmax_sd/(sqrt(n)),
  Jmax_se = Jmax_sd/(sqrt(n))
)


ggplot(aci_parameters_event_summary, aes(x=reorder(event_short,Vcmax)))+
  geom_point(aes(y=Vcmax))+
  geom_errorbar(aes(ymin=Vcmax-(2*Vcmax_se),ymax=Vcmax+(2*Vcmax_se)))

ggplot(aci_parameters_event_summary, aes(x=reorder(event_short,Jmax)))+
  geom_point(aes(y=Jmax))+
  geom_errorbar(aes(ymin=Jmax-(2*Jmax_se),ymax=Jmax+(2*Jmax_se)))

Vcmax_mean <- mean(aci_parameters_list$Vcmax)
Vcmax_sd <- sd(aci_parameters_list$Vcmax)
Vcmax_CV = Vcmax_sd/Vcmax_mean

Jmax_mean <- mean(aci_parameters_list$Jmax)
Jmax_sd <- sd(aci_parameters_list$Jmax)
Jmax_CV = Jmax_sd/Jmax_mean

str(all_trees_wide)


write.csv(aci_parameters_list2, file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/aci_parameters_list.csv")
