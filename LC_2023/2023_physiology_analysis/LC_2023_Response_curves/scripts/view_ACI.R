#script for viewing aci data

library(plantecophys)
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(plantecophys)
library(mgcv)

#read in compiled dataset as a csv

ACI_curves <- read.csv(file = "")


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

