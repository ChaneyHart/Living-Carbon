##Leaf trait analysis 2022

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)
library(emmeans)
library(nlme)
library(lme4)


leaf_area_2022 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/2022_leaf_morphology/LB_Leafarea_2022_fixed.csv", skip =2)
str(leaf_area_2022)
leaf_area_2022 <- na.omit(leaf_area_2022)

growth <- read.csv(file = 'LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv')



leaf_mass_2022 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/2022_leaf_morphology/Leaf_mass_LC_2022.csv")
str(leaf_mass_2022)
leaf_mass_2022$Leaf.3 <- as.numeric(leaf_mass_2022$Leaf.3)
leaf_mass_2022$avg_mass <- as.numeric(leaf_mass_2022$avg_mass)
leaf_mass_2022 <- na.omit(leaf_mass_2022)


SLA_2022 <- inner_join(leaf_area_2022,leaf_mass_2022)
SLA_2022$SLA <- SLA_2022$average_LA/SLA_2022$avg_mass
SLA_2022$LMA <- 1/SLA_2022$SLA

write.csv(SLA_2022, file = "LC_2022/2022_leaf_trait_analysis/2022_leaf_morphology/SLA_2022_raw.csv")




