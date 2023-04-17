#A-Ci data analysis from July 2021

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)



growth_2022 <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
Aci_parameters <- read.csv("2022_physiology_analysis/CO2_response/CO2_response_modeling/Aci_parameters_list.csv")
Aci_parameters$ID <- Aci_parameters$tree

aci_obs_df <- read.csv("2022_physiology_analysis/CO2_response/compiled/ACI_compilation_w_weather.csv")
aci_obs_df$ID <- aci_obs_df$tree


#combine data sets
Aci_growth <- inner_join(Aci_parameters, growth_2022, by = "ID")


aci_obs_subset <- subset(aci_obs_df, select = c(ID,obs,date,hhmmss,instrument,row,column,spad,operator,temp,event,event_short,block,leaf,tree,E,Emm,A,Csetpoint,Ca,Ci,gsw,gbw,gtw,TleafEB,RHcham,VPcham,VPDleaf,PhiPS2,ETR,PhiCO2,Qin,Tleaf))
aci_obs_subset <- left_join(aci_obs_subset, growth_2022, by = "ID")
aci_obs_subset_flr <- subset(aci_obs_subset, PhiPS2 > 0)

###summary graphs

aci_obs_subset_II <- subset(aci_obs_subset, Ci > 0)

ACI_curves_raw <- ggplot(aci_obs_subset_II, aes(x=Ci, y=A, color = event_short.x))+
  geom_point()+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s)")

ggsave('2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/ACI_curves_raw.png', plot = ACI_curves_raw, dpi = 300)

ggplot(aci_obs_subset_II, aes(x=Ci, y=A, color = event_short.x))+
  geom_smooth(span =0.7)+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s)")


ggplot(aci_obs_subset_II, aes(x=Ci, y=A, color = construct2))+
  geom_point()+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s")

ggplot(aci_obs_subset_II, aes(x=Ci, y=A, color = construct2))+
  geom_smooth(span= 0.5)+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s")


#Looking at underlying changes in weather

Air_timecourse <- ggplot(Aci_growth,aes(x=X.x, y = Air_temp))+
  geom_point()+
  geom_path()+
  ylab('Ambient temp (˚C)')+
  xlab('tree')
Air_timecourse
ggsave("2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/Temp_timecourse.png",plot = Air_timecourse, dpi = 300)

VPD_timecourse <- ggplot(Aci_growth,aes(x=X.x, y = Air_VPD))+
  geom_point()+
  geom_path()+
  ylab('Ambient VPD (KPa)')+
  xlab('tree')

ggsave('2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/VPD_timecourse.png', plot = VPD_timecourse, dpi = 300)

#explore covariates

SPAD_covariate <- ggplot(Aci_growth, aes(x = spad, y = A_A_410))+
  geom_point(aes(color = event_short),size = 3)+
  ylab("Ambient assimilation (mol/m^2*s)")+
  xlab("Chollorphyll density (SPAD)")
  
SPAD_covariate
ggsave('2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/SPAD_covariate.png', plot = SPAD_covariate, dpi = 300)

ggplot(Aci_growth, aes(x = Jmax, y = spad))+
  geom_point()

ggplot(Aci_growth, aes(x = J.Vc, y = spad))+
  geom_point()

Air_temp_covariate <- ggplot(Aci_growth,aes(x = Air_temp, y = A_A_410))+
  geom_point(aes(color = event_short),size = 3)+
  ylab("Ambient assimilation (mol/m^2*s)")+
  xlab("Ambient air temp (˚C)")
Air_temp_covariate 

ggsave("2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/Air_temp_covariate.png", plot = Air_temp_covariate, dpi = 300)

Air_VPD_covariate <- ggplot(Aci_growth, aes(x = Air_VPD, y = A_A_410))+
  geom_point(aes(color = event_short),size = 3)+
  ylab("Ambient assimilation (mol/m^2*s)")+
  xlab("Ambient air VPD (KPa)")
Air_VPD_covariate
ggsave("2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/Air_VPD_covariate.png", plot = Air_VPD_covariate, dpi = 300)

ggplot(Aci_growth, aes(x=gsw_A_410, y = Vcmax))+
  geom_point()

ggplot(Aci_growth, aes(x=gsw_A_410, y=A_A_410))+
  geom_point()

A_gsw <- lm(A_A_410 ~ gsw_A_410, Aci_growth)
summary(A_gsw)

#explopratory models

model_test <- lm(Vcmax ~ block + Air_temp + Air_VPD + spad, Aci_growth)
summary(model_test)

#Indicates that environmental factors (air - VPD) and developmental/productivity factors (spad) are significant.
# may want to include date.

model_test2 <- lm(Vcmax ~ block + Air_temp + Air_VPD, Aci_growth)
summary(model_test2)

model_test3 <- lm(Vcmax ~ block + Air_VPD + tempset, Aci_growth)
summary(model_test3)

#need to do  more model selection but will go with #2 for now.
#now analyzing differences between events and constructs

Field_by_event <- group_by(Aci_growth, event_short)
Field_sum <- Field_by_event %>% summarise(
  VI_sd = sd(V497),
  VI_n = dplyr::n(), 
  VI_se = VI_sd / sqrt(VI_n),
  SPAD_sep = mean(Aug_SPAD),
  SPAD_sep_sd = sd(Aug_SPAD),
  SPAD_sep_n = dplyr::n(), 
  SPAD_sep_se = SPAD_sep_sd / sqrt(SPAD_sep_n),
  VcMax = mean(Vcmax),
  VcMax_sd = sd(Vcmax),
  VcMax_n = dplyr::n(), 
  VcMax_se = VcMax_sd / sqrt(VcMax_n),
  Jmax_sd = sd(Jmax),
  Jmax = mean(Jmax),
  Jmax_n = dplyr::n(), 
  Jmax_se = Jmax_sd / sqrt(Jmax_n),
  J_Vc_sd = sd(J.Vc),
  J_Vc = mean(J.Vc),
  J_Vc_n = dplyr::n(),
  J_Vc_se = J_Vc_sd/J_Vc_n
  
)

Vcmax_comp <- ggplot(Aci_growth, aes(x = reorder(event_short, Vcmax), y = Vcmax, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Vcmax (µmol/m^2*s)")

#block effect?
Vcmax_comp2 <- ggplot(Aci_growth, aes(x = reorder(event_short,Vcmax), y = Vcmax, fill = block))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Vcmax (µmol/mol s)")

Vcmax_comp
Vcmax_comp2
#look for outliers 

aov_Vcmax <- aov(Vcmax ~ event_short + block, data = Aci_growth)

plot(aov_Vcmax,which = 1)
plot(aov_Vcmax, which = 2)
#looks good

ggsave("2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/Vcmax_comp.png",plot = Vcmax_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("Vcmax_comp_blockeff.png",plot = Vcmax_comp2, width = 8, height = 5, units = "in", dpi = 300)

##stats######


#####1-way ANOVA w/ #####
anova_Vcmax <- lm(Vcmax ~ event_short + block, data = Aci_growth)
plot(aov_Vcmax, which = 2)
#normality and homogeneity of variance for groups look ok
anova(anova_Vcmax)


#Tukey contrasts####
TukeyHSD(aov_Vcmax)

##Dunnett's comparison w/ escape as control

summary(Aci_growth$construct2)
Aci_growth$construct2 <- as.factor(Aci_growth$construct2)

Vcmax_Dunnets_aov <- aov(Vcmax ~ construct2 + block, data = Aci_growth)
Vcmax_glht <- glht(Vcmax_Dunnets_aov,linfct = mcp(construct2 = "Dunnett"))
summary(Vcmax_glht)


#######mixed model##########
library(car)
library(emmeans)
library(nlme)

mmodel_0 <- lm(Vcmax ~ block + Air_temp + Air_VPD + spad, Aci_growth)
summary(mmodel_0)
vif(mmodel_0)
#air temp and Air VPD collinear
# remove VPD from model, air temp sufficient

mmodel_1 <- lm(Vcmax ~ block + Air_temp + spad, Aci_growth)
summary(mmodel_1)
vif(mmodel_1)

Vcmax_mod1 <- lme(Vcmax ~ block + construct2, random = ~1|Air_temp, data = Aci_growth)
summary(Vcmax_mod1)

plot(fitted(Vcmax_mod1), residuals(Vcmax_mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(Vcmax_mod1)); qqline(residuals(Vcmax_mod1))
#residuals NOT randomly distributed, fishy

##############################3
#Jmax

Jmax_comp <- ggplot(Aci_growth, aes(x = reorder(event_short, Jmax), y = Jmax, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Jmax (µmol/mol s)")

#block effect?
Jmax_comp2 <- ggplot(Aci_growth, aes(x = reorder(event_short,Jmax), y = Jmax, fill = block))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Jmax (µmol/mol s)")


Jmax_comp
Jmax_comp2
#look for outliers 

aov_Jmax <- aov(Jmax ~ event_short + block, data = Aci_growth)

plot(aov_Jmax,which = 1)
#looks good

ggsave("2022_physiology_analysis/CO2_response/CO2_response_graph_analysis/Jmax_comp.png",plot = Jmax_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("Jmax_comp_blockeff.png",plot = Jmax_comp2, width = 8, height = 5, units = "in", dpi = 300)

##stats######


#####1-way ANOVA w/ #####
anova_Jmax <- lm(Jmax ~ event_short + block, data = Aci_growth)
anova(anova_Jmax)
#There is moderate evidence of a difference that there is a difference in mean height for event (p-value: 0.09066)

#Tukey contrasts####
TukeyHSD(aov_Jmax)

##Dunnett's comparison w/ escape as control

summary(Aci_growth$construct2)
Aci_growth$construct2 <- as.factor(Aci_growth$construct2)

Jmax_Dunnets_aov <- aov(Jmax ~ construct2 + block, data = Aci_growth)
Jmax_glht <- glht(Jmax_Dunnets_aov,linfct = mcp(construct2 = "Dunnett"))
summary(Jmax_glht)


#########
#ratio
J.Vc_comp <- ggplot(Aci_growth, aes(x = reorder(event_short, J.Vc), y = J.Vc, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Jmax/Vcmax ")
J.Vc_comp

#block effect?
J.Vc_comp2 <- ggplot(Aci_growth, aes(x = reorder(event_short,J.Vc), y = J.Vc, fill = block))+
  geom_boxplot()+
  xlab("Event")+
  ylab("J/Vc")


J.Vc_comp
J.Vc_comp2
#look for outliers 

aov_J.Vc <- aov(J.Vc ~ event_short + block, data = Aci_growth)

plot(aov_J.Vc,which = 1)
#looks good

ggsave("J.Vc_comp.png",plot = J.Vc_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("J.Vc_comp_blockeff.png",plot = J.Vc_comp2, width = 8, height = 5, units = "in", dpi = 300)

##stats######


#####1-way ANOVA w/ #####
anova_J.Vc <- lm(J.Vc ~ event_short + block, data = Aci_growth)
anova(anova_J.Vc)
#There is moderate evidence of a difference that there is a difference in mean height for event (p-value: 0.09066)

#Tukey contrasts####
TukeyHSD(aov_J.Vc)

##Dunnett's comparison w/ escape as control

summary(growth_phen_compare_inner_Dunnetts$construct)

J.Vc_Dunnets_aov <- aov(H356 ~ construct + block, data = growth_phen_compare_inner_Dunnetts)
J.Vc_glht <- glht(J.Vc_Dunnets_aov,linfct = mcp(construct = "Dunnett"))
summary(J.Vc_glht)

#
#correlating growth to phy
growth_phys$V419 <- as.numeric(growth_phys$V419)
growth_phys$Vcmax <- as.numeric(growth_phys$Vcmax)

growth_phys <- ggplot(Aci_growth,aes(x = V419, y =Vcmax), fill = event_short)+
  geom_point()+
  xlab('volume index (mm3)')+
  ylab('Vcmax (µmol/mol s)')

growth_phys_model <- lm(Vcmax~H419, Aci_growth)
summary(growth_phys_model)
# r2 near 0
growth_phys

growth_phys2 <- ggplot(Aci_growth, aes(x = H419, y = Vcmax), fill = event_short)+
  geom_point()+
  xlab('height(mm)')+
  ylab('Vcmax (µmol/mol s)')

growth_phys_3 <- ggplot(Aci_growth, aes(x = H419, y = Jmax), fill = event_short)+
  geom_point()

growth_phys2                      
growth_phys_3
