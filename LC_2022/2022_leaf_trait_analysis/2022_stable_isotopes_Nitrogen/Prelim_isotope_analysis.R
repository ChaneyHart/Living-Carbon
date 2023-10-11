#Prelim isotope analysis
#5-28-22
#Explore relationship between %N and SPAD

library(dplyr)
library(scales)
library(tidyr)
library(ggplot2)



raw_isotope <- read.csv("2022_leaf_trait_analysis/2022_leaf_morphology/Prelim_isotope_R.csv")

#clean up
raw_isotope
raw_isotope <- rename(raw_isotope, "ID" = "Sample")
raw_isotope$GreenessScore <- as.factor(raw_isotope$GreenessScore)

quickplot(Event, wt.N, data = raw_isotope)
quickplot(Event, Spad_Sep, data = raw_isotope,color = GreenessScore)
quickplot(wt.N,Spad_Sep, data = raw_isotope, color = Event)
quickplot(wt.N, d13C, data = raw_isotope, color = Block)
quickplot(d13C,Spad_Sep, data = raw_isotope, color = Event)

quickplot(wt.N, Spad_sep, data= raw_isotope, color = GreenessScore)

#SPAD and %N
ggplot(data = raw_isotope, aes(x = wt.N, y = Spad_Sep, color = GreenessScore))+
  geom_point()
  
growth_summ <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")  
growth_summ <- rename(growth_summ,"ID" = "id")

isotope_plus_growth <- inner_join(raw_isotope, growth_summ, by = "ID")
sapply(isotope_plus_growth,mode)
isotope_plus_growth$Spad_Sep <- as.numeric(isotope_plus_growth[,13])
isotope_plus_growth$GreenessScore <- as.factor(isotope_plus_growth$GreenessScore)

pairs(isotope_plus_growth[,c(6:8)])
pairs(isotope_plus_growth[,c(12,13,25,26)])

SPAD_N_anova <- lm(Spad_Sep ~ wt.N, data = isotope_plus_growth)
anova(SPAD_N_anova)

#iWUE
prelim_iWUE <- ggplot(data = isotope_plus_growth, aes(x = reorder(Event,iWUE), y = iWUE, fill = construct))+
  geom_boxplot()+
  labs(y = "iWUE (µmol CO2/mmol H20)", x = "Event")

prelim_iWUE

##stats###

iWUEanova <- lm(iWUE ~ Event + Block, data = isotope_plus_growth)
anova(iWUEanova)

library(emmeans)

emmeans(iWUEanova, specs = pairwise ~ Event)
mean(isotope_plus_growth$iWUE)

ggsave(filename = "prelim_iWUE.png", plot = prelim_iWUE, dpi = 300)

WUE_comp <- ggplot(data = isotope_growth, aes(x = reorder(Event,iWUE), y = iWUE, fill = construct))+
  geom_boxplot()+
  labs(x = "Event", y = "iWUE (µmol CO2/mmol H20)")
WUE_comp

WUEcomp2 <- ggplot(data = isotope_growth, aes(x = reorder(construct,iWUE), y = iWUE, fill = block))+
  geom_boxplot()+
  labs(x = "Event", y = "Ci (µmol CO2/mmol H20)")
WUEcomp2

ggsave("WUE_event_comp.png", plot = WUE_comp, width = 6, height = 3, dpi = 300, units = "in")

iWUEanova <- lm(iWUE ~ Event + Block, data = isotope_growth)
anova(iWUEanova)

growth_0722 <- read.csv("LC_June2022/DBH_height_timeline_exluded.csv")
growth_0722$R298 <- as.numeric(growth_0722$R298)
growth_0722$R334 <- as.numeric(growth_0722$R334)
growth_0722$R356 <- as.numeric(growth_0722$R356)

isotope_growth <- inner_join(raw_isotope, growth_0722, by = "ID")
isotope_growth$Event <- as.factor(isotope_growth$Event)


iWUE_VI <- ggplot(data = isotope_growth, aes(x= iWUE, y = V144, color = construct, shape = event))+
  geom_point()+
  labs(x = "iWUE", y = "December 2021 VI (mm^3)")

iWUE_VI

lmiWUE_VI <- lm(V144 ~ iWUE + Event, data = isotope_growth)
summary(lmiWUE_VI)

ggsave("iWUE_VI.png", plot = iWUE_VI, width = 6, height = 3, dpi = 300, units = "in")
