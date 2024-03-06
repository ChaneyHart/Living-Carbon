library(emmeans)
library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)
#import cleaned AQ data

#import photosynthetic parameters
AQ_parameters <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/AQ_parameters_list.csv")

#stats


# Create tiers factor
AQ_parameters <- AQ_parameters %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "intermediate",
  event_short == "13-15E" | event_short == "2H" ~ "high",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))


AQ_parameters$sample_date <- as.factor(AQ_parameters$sample_date)
# summary stats
mean_k_sat <- mean(AQ_parameters$k_sat)
std_k_sat <- sd(AQ_parameters$k_sat)
CV_k_sat <- (std_k_sat/mean_k_sat)*100

library(lmerTest)

k_sat_mod <- lmer(k_sat ~ Class + (1|sample_date), data = subset(AQ_parameters,Class !="WT"))

AQ_parameters <- subset(AQ_parameters,Class != "WT") %>% mutate(
  k_sat_resid = resid(k_sat_mod),
  k_sat_fitted = fitted(k_sat_mod)
)

#examine residuals
ggplot(AQ_parameters,aes(x=k_sat_fitted,y=k_sat_resid))+
  geom_point()

ggplot(AQ_parameters,aes(x=Class,y=k_sat_resid))+
  geom_point()

ggplot(AQ_parameters,aes(x=Class,y=k_sat_resid))+
  geom_point()

#pairwise comparisons


emmeans_k_sat <- emmeans(k_sat_mod, specs = pairwise ~Class)
emmeans_k_sat

# graph modelled tier means
Class_k_sat_effect <- as.data.frame(emmeans_k_sat$emmeans)
Class_Jmax_effect$modelled_k_sat <- (Class_k_sat_effect$emmean)

my_colors3 <- c("gray0","chartreuse4","indianred3")
modelled_k_sat_plot <- ggplot(Class_k_sat_effect, aes(x=Class,y=emmean,fill=Class))+
  geom_bar(stat = "identity")+
  ylab("A sat (µmol/m2/sec)")+
  geom_errorbar(aes(ymin=emmean-2*SE,ymax=emmean+2*SE),linewidth=0.5,width=0.5,color="gray60")+
  scale_fill_manual(values = my_colors3)+
  theme_bw()

modelled_k_sat_plot


####
phi_J_mod <- lmer(phi_J ~ Class + (1|sample_date), data = subset(AQ_parameters,Class !="WT"))

AQ_parameters <- subset(AQ_parameters,Class != "WT") %>% mutate(
  phi_J_resid = resid(phi_J_mod),
  phi_J_fitted = fitted(phi_J_mod)
)

#examine residuals
ggplot(AQ_parameters,aes(x=phi_J_fitted,y=phi_J_resid))+
  geom_point()

#maybe some non-constant variance

ggplot(AQ_parameters,aes(x=Class,y=phi_J_resid))+
  geom_point()

ggplot(AQ_parameters,aes(x=sample_date,y=phi_J_resid))+
  geom_point()

#pairwise comparisons


emmeans_phi_J <- emmeans(phi_J_mod, specs = pairwise ~Class)
emmeans_phi_J

# graph modelled tier means
Class_phi_J_effect <- as.data.frame(emmeans_phi_J$emmeans)
Class_Jmax_effect$modelled_phi_J <- (Class_phi_J_effect$emmean)

my_colors3 <- c("gray0","chartreuse4","indianred3")
modelled_phi_J_plot <- ggplot(Class_phi_J_effect, aes(x=Class,y=emmean,fill=Class))+
  geom_bar(stat = "identity")+
  ylab("Quantum efficiency (mol CO2/ mol J)")+
  geom_errorbar(aes(ymin=emmean-2*SE,ymax=emmean+2*SE),linewidth=0.5,width=0.5,color="gray60")+
  scale_fill_manual(values = my_colors3)+
  theme_bw()

modelled_phi_J_plot

