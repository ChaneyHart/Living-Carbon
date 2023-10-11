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


leaf_area_2022 <- read.csv(file = "2022_leaf_trait_analysis/2022_leaf_morphology/LB_Leafarea_2022_fixed.csv", skip =2)

leaf_area_2022 <- na.omit(leaf_area_2022)

growth <- read.csv(file = '2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv')

leaf_area_dat <- inner_join(leaf_area_2022, growth, by = "ID")


##plot
ggplot(leaf_area_dat, aes(x=event_short, y=average_LA, fill = construct2))+
  geom_boxplot()
