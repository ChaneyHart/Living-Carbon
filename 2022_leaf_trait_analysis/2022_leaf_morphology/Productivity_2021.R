#revisiting %N data

library(dplyr)
library(scales)
library(tidyr)
library(ggplot2)
library(corrplot)
library(emmeans)
library(nlme)

raw_isotope <- read.csv("2022_leaf_trait_analysis/2022_leaf_morphology/Prelim_isotope_R.csv")
raw_isotope <- dplyr::rename(raw_isotope, "ID" = "Sample")
raw_isotope$GreenessScore <- as.factor(raw_isotope$GreenessScore)
growth_summ <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv") 
SLA <- read.csv("2021_LMA.csv")
SLA <- subset(SLA, select = c(ID, SLA.cm2.g.))


Nitrogen_dat <- inner_join(raw_isotope,growth_summ, by ="ID")
Nitrogen_dat <- inner_join(Nitrogen_dat, SLA, by="ID")
Nitrogen_dat <- subset(Nitrogen_dat, SPAD_sep2021 > 0)
Nitrogen_dat$MLA_2021 <- as.numeric(Nitrogen_dat$MLA_2021)
nitrogen_pairs <- subset(Nitrogen_dat, select = c(wt.N,dec_VI,SLA.cm2.g., SPAD_sep2021))
nitrogen_pairs <- nitrogen_pairs %>% dplyr::rename(percent_N = wt.N, logVI = dec_VI, SLA = SLA.cm2.g., Chl_density = SPAD_sep2021)
nitrogen_pairs <- subset(nitrogen_pairs, Chl_density > 0)
nitrogen_pairs <- subset(nitrogen_pairs, SLA > 0)


nitrogen_pairs$SLA <- log(nitrogen_pairs$SLA)
nitrogen_pairs$logVI <- log(nitrogen_pairs$logVI)
nitrogen_pairs$precent_N <- log(nitrogen_pairs$precent_N)
nitrogen_pairs$SPAD <- log(nitrogen_pairs$SPAD)

#correlation plots
install.packages("GGally")
library(GGally)
str(nitrogen_pairs)
prod_plot <- ggpairs(nitrogen_pairs)
ggsave(filename = "prod_plot.png", plot = prod_plot, dpi = 300)
cor <- cor(nitrogen_pairs)
corrplot(cor, method = "number")

#comparing constructs
ggplot(Nitrogen_dat, aes(x=Event, y = wt.N, fill = construct2))+
  geom_boxplot()+
  xlab("Event")+
  ylab("%Nitrogen (g/g)")
  
ggplot(Nitrogen_dat, aes(x=construct2, y = wt.N, fill = construct2))+
  geom_boxplot()+
  xlab("Event")+
  ylab("%Nitrogen (g/g)")
  
N_anova <- lm(wt.N ~ Event.x, Nitrogen_dat)
summary(N_anova)
anova(N_anova)

N_anova2 <- lm(wt.N ~ construct2, Nitrogen_dat)
summary(N_anova2)
anova(N_anova2)

ggplot(Nitrogen_dat, aes(x=Event.x, y=SLA.cm2.g., fill = construct2))+
  geom_boxplot()

ggplot(Nitrogen_dat, aes(x=Event.x, y=SPAD_sep2021, fill=construct2))+
  geom_boxplot()

#is the relationship between %N and SPAD significant
N_spad_model <- lm(wt.N ~ SPAD_sep2021, Nitrogen_dat)
summary(N_spad_model)

N_SPAD <- ggplot(Nitrogen_dat, aes(x=SPAD_sep2021, y = wt.N))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("SPAD (september 2021)")+
  ylab("% Nitrogen (g/g)")

ggsave(filename = "N_SPAD.png", plot = N_SPAD, dpi = 300)

N_spad_model_par <- lm(wt.N ~ SPAD_sep2021 + construct2, Nitrogen_dat)
summary(N_spad_model_par)

N_spad_model_int <- lm(wt.N ~ SPAD_sep2021 + construct2 + SPAD_sep2021*construct2, Nitrogen_dat)
summary(N_spad_model_int)
