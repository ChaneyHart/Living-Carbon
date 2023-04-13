#A-Ci data analysis from July 2021

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)



growth_0912 <- read.csv("2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
Aci_parameters <- read.csv("2022_physiology_analysis/CO2_response/CO2_response_modeling/Aci_parameters_list.csv")

aci_obs_df <- read.csv("2022_physiology_analysis/CO2_response/compiled/ACI_compilation_w_weather.csv")
aci_obs_df$ID <- aci_obs_df$tree

#combine data sets
Aci_growth <- inner_join(Aci_parameters, growth_0912, by = "ID")


aci_obs_subset <- subset(aci_obs_df, select = c(ID,obs,date,hhmmss,instrument,row,column,spad,operator,temp,event,event_short,block,leaf,tree,E,Emm,A,Csetpoint,Ca,Ci,gsw,gbw,gtw,TleafEB,RHcham,VPcham,VPDleaf,PhiPS2,ETR,PhiCO2,Qin,Tleaf))
aci_obs_subset <- left_join(acid_obs_subset, growth_0912, by = "ID")
aci_obs_subset_flr <- subset(aci_obs_subset, PhiPS2 > 0)

###summary graphs

all_trees_df_subset <- subset(all_trees_df_subset, Ci > 0)

ggplot(all_trees_df_subset, aes(x=Ci, y=A, color = event_short.x))+
  geom_point()+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s")

ggplot(all_trees_df_subset, aes(x=Ci, y=A, color = event_short.x))+
  geom_smooth(span =0.7)+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s)")



ggplot(all_trees_df_subset, aes(x=Ci, y=A, color = construct2))+
  geom_point()+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s")

ggplot(all_trees_df_subset, aes(x=Ci, y=A, color = construct2))+
  geom_smooth(span= 0.5)+
  xlab("Ci (ppm)")+
  ylab("A (µmol/m^s*s")


Aci_growth <- inner_join(Aci_parameters, growth_0912, by = "ID")

ggplot(Aci_growth,aes(x=obs, y = Air_temp))+
  geom_point()+
  geom_path()+
  ylab('Ambient temp (˚C)')+
  xlab('tree')

ggplot(Aci_growth,aes(x=obs, y = Air_VPD))+
  geom_point()+
  geom_path()+
  ylab('Ambient VPD (KPa)')+
  xlab('tree')


#explore the data

ggplot(Aci_growth, aes(x = SPAD, y = A_A_410))+
  geom_point()+
  ylab("ambient assimilation (mol/m^2*s)")
  

ggplot(Aci_growth, aes(x = Jmax, y = SPAD))+
  geom_point()

ggplot(Aci_growth, aes(x = J.Vc, y = SPAD))+
  geom_point()

ggplot(Aci_growth,aes(x = Air_temp, y = A_A_410))+
  geom_point()+
  ylab("ambient assimilation (mol/m^2*s)")+
  xlab("ambient air temp (˚C)")

ggplot(Aci_growth, aes(x = Air_VPD, y = A_A_410))+
  geom_point()+
  ylab("ambient assimilation (mol/m^2*s)")+
  xlab("ambient air VPD (KPa)")

ggplot(Aci_growth, aes(x=gsw_A_410, y = Vcmax))+
  geom_point()

ggplot(Aci_growth, aes(x=gsw_A_410, y=A_A_410))+
  geom_point()

A_gsw <- lm(A_A_410 ~ gsw_A_410, Aci_growth)
summary(A_gsw)

#explopratory model

model_test <- lm(Vcmax ~ block + Air_temp + Air_VPD + SPAD + tempset, Aci_growth)
summary(model_test)

model_test2 <- lm(Vcmax ~ block + Air_temp + Air_VPD, Aci_growth)
summary(model_test2)

model_test3 <- lm(Vcmax ~ block + Air_VPD + tempset, Aci_growth)
summary(model_test3)

#need to do  more model selection but will go with #2 for now.
#now analyzing differences between events and constructs

Field_sum <- Field_by_event %>% summarise(
  VI_sd = sd(V419),
  VI_n = n(), 
  VI_se = VI_sd / sqrt(VI_n),
  SPAD_sep = mean(Aug_SPAD),
  SPAD_sep_sd = sd(Aug_SPAD),
  SPAD_sep_n = n(), 
  SPAD_sep_se = SPAD_sep_sd / sqrt(SPAD_sep_n),
  VcMax = mean(Vcmax),
  VcMax_sd = sd(Vcmax),
  VcMax_n = n(), 
  VcMax_se = VcMax_sd / sqrt(VcMax_n),
  Jmax_sd = sd(Jmax),
  Jmax = mean(Jmax),
  Jmax_n = n(), 
  Jmax_se = Jmax_sd / sqrt(Jmax_n),
  J_Vc_sd = sd(J.Vc),
  J_Vc = mean(J.Vc),
  J_Vc_n = n(),
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
#looks good

ggsave("Vcmax_comp.png",plot = Vcmax_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("Vcmax_comp_blockeff.png",plot = Vcmax_comp2, width = 8, height = 5, units = "in", dpi = 300)

##stats######


#####1-way ANOVA w/ #####
anova_Vcmax <- lm(Vcmax ~ event_short + block, data = Aci_growth)
anova(anova_Vcmax)
#There is moderate evidence of a difference that there is a difference in mean height for event (p-value: 0.09066)

#Tukey contrasts####
TukeyHSD(aov_Vcmax)

##Dunnett's comparison w/ escape as control

summary(Aci_growth$construct2)
Aci_growth$construct2 <- as.factor(Aci_growth$construct2)

Vcmax_Dunnets_aov <- aov(Vcmax ~ construct2 + block, data = Aci_growth)
Vcmax_glht <- glht(Vcmax_Dunnets_aov,linfct = mcp(construct2 = "Dunnett"))
summary(Vcmax_glht)

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

ggsave("Jmax_comp.png",plot = Jmax_comp, width = 8, height = 5, units = "in", dpi = 300)
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
