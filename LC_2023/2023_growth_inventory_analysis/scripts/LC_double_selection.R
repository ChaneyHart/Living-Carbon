#Living Carbon Double selection

library(tidyverse)
#install.packages("caret")
library(caret)
#install.packages("leaps")
library(leaps)
library(MASS)
library(lmerTest)

#read in dat 

full_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/processed_growth_dat.csv")

#first step, regress Y (biomass) on all covariates

covariates <- subset(full_dat, select = -c(X,ID,row,column,event_short,block,Class,H801,D801,D801,geometry))

#perform model selection of all of the covariates on Y (V801)



cov_lm <- lm(V801 ~., data=covariates)

cov_step <- stepAIC(cov_lm, direction = "both")

summary(cov_step)

cov_model <- lm(V801 ~ D49 + V49 + SPAD_June_2022 + SPAD_sep_2022 + 
                  SPAD_Aug_2023 + SPAD_Sep_2023 + nir_mean + re_mean + MGRVI_mean + 
                  NDVI_mean + NDRE_mean + SAVI_mean + OSAVI_mean + shadow_area + 
                  branchiness + area_avg, data = covariates)

#examine collinearity




