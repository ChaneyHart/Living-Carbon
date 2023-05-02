#LC growth analysis update - March 24 2022
##This script was used to perform analysis of the Living Carbon trial growth data up to this point in March 2022
#results in LC_update_lab_meeting_9_12


library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)



#bring in growth data w/ info about each tree
# Data has been cleaned in 2022_field_season_data_cleaning.R
# Negative growth and missing data was imputted for diameter and height. No outliers were removed
# CT1 was excluded as were trees that were damaged and replaced or died.

#read in data
growth <- read.csv("2022_growth&inventory_analylsis/growth_data_cleaning/LC_3_22_growth_data_cleaned.csv")
#check data types
str(growth)
growth$D385 <- as.numeric(growth$D385)

#calculate volume

growth$V49 = ((growth$D49^2)*growth$H49)
growth$V144 = ((growth$D144^2)*growth$H144)
growth$V299 = ((growth$D299^2)*growth$H299)
growth$V335 = ((growth$D335^2)*growth$H335)
growth$V357 = ((growth$D357^2)*growth$H357)
growth$V385 = ((growth$D385^2)*growth$H385)
growth$V421 = ((growth$D421^2)*growth$H421)
growth$V497 = ((growth$D497^2)*growth$H497)

growth$V49_cm = ((growth$D49^2)*growth$H49)/1000
growth$V144_cm = ((growth$D144^2)*growth$H144)/1000
growth$V299_cm = ((growth$D299^2)*growth$H299)/1000
growth$V335_cm = ((growth$D335^2)*growth$H335)/1000
growth$V357_cm = ((growth$D357^2)*growth$H357)/1000
growth$V385_cm = ((growth$D385^2)*growth$H385)/1000
growth$V421_cm = ((growth$D421^2)*growth$H421)/1000
growth$V497_cm = ((growth$D497^2)*growth$H497)/1000


#Convert dataframe to long format to look at changes over time for each tree
Height_long <- pivot_longer(growth, cols = c(H49,H144,H299,H335,H357,H385,H421,H497), names_to = "Days", values_to = "Height")
Height_long <- Height_long %>% mutate(Days = as.numeric(substr(Height_long$Days,2,4)))
Height_long$meters <- Height_long$Height/1000

Diam_long <- pivot_longer(growth, cols = c(D49,D144,D299,D335,D357,D385,D421,D497), names_to = "Days", values_to = "Diameter")
Diam_long <- Diam_long %>% mutate(Days = as.numeric(substr(Diam_long$Days,2,4)))

Volume_long <- pivot_longer(growth, cols = c(V49,V144,V299,V335,V357,V385,V421,V497), names_to = "Days", values_to = "Volume")
Volume_long <- Volume_long %>% mutate(Days = as.numeric(substr(Volume_long$Days,2,4)))
Volume_long$cm3 <- Volume_long$Volume/1000


ht_full <- ggplot(Height_long, aes(x = Days, y = meters))+
  geom_point(aes(color = event_short), show.legend = FALSE)+
  xlab("Days since planting")+
  ylab("Height (m)")


ht_full
ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline.png", plot = ht_full, width = 6, height = 3, units = "in", dpi = 300)

d_full <- ggplot(Diam_long, aes(x = Days, y = Diameter))+
  geom_point(aes(color = event_short), show.legend = FALSE)+
  xlab("days since planting")+
  ylab("Diameter (mm)")


d_full
ggsave("2022_growth&inventory_analylsis/growth_graph/d_full.png", plot = d_full, width = 6, height = 3, units = "in", dpi = 300)


vol_full <- ggplot(Volume_long, aes(x = Days, y = (cm3)))+
  geom_point(aes(color = event_short), show.legend = FALSE)+
  xlab("days since planting")+
  ylab("Volume (cubic cm)")



vol_full
ggsave("2022_growth&inventory_analylsis/growth_graph/vol_full.png", plot = vol_full, width = 6, height = 3, units = "in", dpi = 300)

### same thing but broken down by event means######


#Height###

Height_long_II <- Height_long %>% group_by(event_short,Days,block)

Height_long_II_summary <- Height_long_II %>% dplyr::summarise(
  height = mean(meters),
  n = n(),
  height_sd = sd(meters),
  height_se = height_sd/(sqrt(n))
)

Height_long_III <- Height_long %>% group_by(event_short,Days)

Height_long_III_summary <- Height_long_III %>% dplyr::summarise(
  height = mean(meters),
  n = n(),
  height_sd = sd(meters),
  height_se = height_sd/(sqrt(n))
)

ht_full_event <- ggplot(Height_long_III_summary, aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE)+
  xlab("Days since planting")+
  ylab("Height (m)")
#large block
ht_full_event_lb <- ggplot(subset(Height_long_II_summary, block == "large"), aes(x = Days, y = height))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("height (m)")
##small block
ht_full_event_sb <- ggplot(subset(Height_long_II_summary, block == "small"), aes(x = Days, y = height))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("height (m)")

  

ht_full_event
ht_full_event_lb
ht_full_event_sb
ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event.png", plot = ht_full_event, width = 6, height = 4, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_lb.png", plot = ht_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_sb.png", plot = ht_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)



###diameter by event #####

diam_long_II <- Diam_long %>% group_by(event_short,Days, block)

diam_long_II_summary <- diam_long_II %>% dplyr::summarise(
  diam = mean(Diameter),
  n = n(),
  diameter_sd = sd(Diameter),
  diameter_se = diameter_sd/(sqrt(n))
)

diam_long_III <- Diam_long %>% group_by(event_short,Days)

diam_long_III_summary <- diam_long_III %>% dplyr::summarise(
  diam = mean(Diameter),
  n = n(),
  diameter_sd = sd(Diameter),
  diameter_se = diameter_sd/(sqrt(n))
)


#all
Diam_full_event <- ggplot(diam_long_III_summary, aes(x = Days, y = diam))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

#large block
Diam_full_event_lb <- ggplot(subset(diam_long_II_summary, block == "large"), aes(x = Days, y = diam))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("diameter (mm)")
##small block
Diam_full_event_sb <- ggplot(subset(diam_long_II_summary, block == "small"), aes(x = Days, y = diam))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

Diam_full_event
Diam_full_event_lb
Diam_full_event_sb
ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event.png", plot = Diam_full_event, width = 6, height = 3, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_lb.png", plot = Diam_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_sb.png", plot = Diam_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)


#####volume by event ###########

volume_long_II <- Volume_long %>% group_by(event_short,Days, block)
volume_long_III <- Volume_long %>% group_by(event_short,Days)

volume_long_II_summary <- volume_long_II %>% dplyr::summarise(
  vol = mean(cm3),
  n = n(),
  vol_sd = sd(cm3),
  vol_se = vol_sd/(sqrt(n))
)

volume_long_III_summary <- volume_long_III %>% dplyr::summarise(
  vol = mean(cm3),
  n = n(),
  vol_sd = sd(cm3),
  vol_se = vol_sd/(sqrt(n))
)


vol_full_event <- ggplot(volume_long_III_summary, aes(x = Days, y = vol))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")

#large block
vol_full_event_lb <- ggplot(subset(volume_long_II_summary, block == "large"), aes(x = Days, y = vol))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")

#small block

vol_full_event_sb <- ggplot(subset(volume_long_II_summary, block == "small"), aes(x = Days, y = vol))+
  geom_jitter(aes(color = event_short), show.legend = TRUE, width = 3)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")


vol_full_event
vol_full_event_lb
vol_full_event_sb
ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event.png", plot = vol_full_event, width = 6, height = 3, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_lb.png", plot = vol_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_sb.png", plot = vol_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)



#Summary statistics
#2021 height growth
(mean(growth$H144) - mean(growth$H49))/mean(growth$H49)
#2021 diameter growth
(mean(growth$D144) - mean(growth$D49))/mean(growth$D49)
#2021 volume growth
(mean(growth$V144) - mean(growth$V49))/mean(growth$V49)


#2022 height growth
mean(growth$H497)
H497_se <- sqrt(var(growth$H497)/(length(growth$H497)))

mean(growth$H299)
H299_se <- sqrt(var(growth$H299)/(length(growth$H299)))
var(growth$H299)
length(growth$H299)


ht_percent_22_upper <- (((mean(growth$H497)+(2*H497_se)) - (mean(growth$H299)+(2*H299_se)))/(mean(growth$H299)+(2*H299_se)))*100
ht_percent_22_lower <- (((mean(growth$H497)-(2*H497_se)) - (mean(growth$H299)-(2*H299_se)))/(mean(growth$H299)-(2*H299_se)))*100
ht_percent_22 <- (((mean(growth$H497)) - (mean(growth$H299)))/(mean(growth$H299)))*100

#2022 diameter growth
mean(growth$D497)
D497_se <- sqrt(var(growth$D497)/(length(growth$D497)))

mean(growth$D299)
D299_se <- sqrt(var(growth$D299)/(length(growth$D299)))
var(growth$D299)
length(growth$D299)


diam_percent_22_upper <- (((mean(growth$D497)+(2*D497_se)) - (mean(growth$D299)+(2*D299_se)))/(mean(growth$D299)+(2*D299_se)))*100
diam_percent_22_lower <- (((mean(growth$D497)-(2*D497_se)) - (mean(growth$D299)-(2*D299_se)))/(mean(growth$D299)-(2*D299_se)))*100
diam_percent_22 <- (((mean(growth$D497)) - (mean(growth$D299)))/(mean(growth$D299)))*100

#2022 volume growth
mean(growth$V497_cm)
V497_se <- sqrt(var(growth$V497_cm)/(length(growth$V497_cm)))
V299_se <- sqrt(var(growth$V299_cm)/(length(growth$V299_cm)))



v_percent_22_upper <- (((mean(growth$V497_cm)+(2*V497_se)) - (mean(growth$V299_cm)+(2*V299_se)))/(mean(growth$V299_cm)+(2*V299_se)))*100
v_percent_22_lower <- (((mean(growth$V497_cm)-(2*V497_se)) - (mean(growth$V299_cm)-(2*V299_se)))/(mean(growth$V299_cm)-(2*V299_se)))*100
v_percent_22 <- (((mean(growth$V497_cm)) - (mean(growth$V299_cm)))/(mean(growth$V299_cm)))*100


(median(growth$V497_cm) - median(growth$V299_cm))/(median(growth$V299_cm))

mean(growth$V299)

###statistical analysis##############
###First go at it by Chaney####
###simple models#####
######Growth analysis##########

###Height and DBH at end of season####




#convert to meters
growth$H497m <- growth$H497/1000

#plot
growth2022_ht_comp <- ggplot(growth, aes(x = reorder(event_short,H497m), y = H497m, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("2022 growing season height(m)")

#block effect?
growth2022_ht_comp2 <- ggplot(growth, aes(x = reorder(construct,H497), y = H497m, fill = block))+
  geom_boxplot()+
  xlab("Event")+
  ylab("2022 growing season height(mm)")


growth2022_ht_comp
growth2022_ht_comp2
#look for outliers 

aov_2022_ht <- aov(H497 ~ event + block, data = growth)
plot(aov_2022_ht,which = 1)
#Identified LCOR-445 as not measured in Sep. Due to beehive?

ggsave("2022_growth&inventory_analylsis/growth_graph/2022_growing_season_ht_comp.png",plot = growth2022_ht_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("2022_growing_season_ht_comp_blockeff.png",plot = growth2022_ht_comp2, width = 8, height = 5, units = "in", dpi = 300)

##stats######


#####1-way ANOVA w/ #####
anova_2022_ht <- lm(H497 ~ event + block, data = growth)
anova(anova_2022_ht)
#There is moderate evidence of a difference that there is a difference in mean height for event (p-value: 0.09066)

#Tukey contrasts####
TukeyHSD(aov_2022_ht)

##Dunnett's comparison w/ escape as control

summary(growth$construct)
growth$construct <- as.factor(growth$construct)

growth_Dunnetts <- growth
growth_Dunnetts$construct <- relevel(growth_Dunnetts$construct, "Escape")
summary(growth_Dunnetts$construct)

growth2022_Dunnets_aov <- aov(H497 ~ construct + block, data = growth_Dunnetts)
growth2022_glht <- glht(growth2022_Dunnets_aov,linfct = mcp(construct = "Dunnett"))
summary(growth2022_glht)


########################################Diameter


growth2022_diam_comp <- ggplot(growth, aes(x = reorder(event_short,D497), y = D497, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("2022 growing season diameter(mm)")

growth2022_diam_comp2 <- ggplot(growth, aes(x = reorder(construct,D497), y = D497, fill = block))+
  geom_boxplot()+
  xlab("Event")+
  ylab("Sep diameter(mm)")

growth2022_diam_comp
growth2022_diam_comp2

#look for outliers

aov_growth2022_diam <- aov(D497 ~ event + block, data = growth)
plot(aov_growth2022_diam,which = 1)
##No missing values


growth2022_diam_comp

ggsave("2022_growth&inventory_analylsis/growth_graph/Growth2022_diam_comp.png",plot = growth2022_diam_comp, width = 8, height = 5, units = "in", dpi = 300)
ggsave("Growth2022_diam_comp_blockeff.png",plot = growth2022_diam_comp2, width = 8, height = 5, units = "in", dpi = 300)

####stats#####
anova_growth2022_diam <- lm(D497 ~ event + block, data = growth)

#event comparison
anova(anova_growth2022_diam)

#######Tukey contrasts#######
TukeyHSD(aov_growth2022_diam)

##Dunnett's comparison w/ escape as control

growth2022_diam_Dunnets_aov <- aov(D497 ~ construct + block, data = growth_Dunnetts)
Sep_diam_glht <- glht(growth2022_diam_Dunnets_aov,linfct = mcp(construct = "Dunnett"))
summary(Sep_diam_glht)

#######
growth2022_vol_comp <- ggplot(growth, aes(x = reorder(event_short,V497_cm), y = V497_cm, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("2022 growing season volume (cubic cm)")

ggsave('2022_growth&inventory_analylsis/growth_graph/growth2022_vol_comp.png', plot = growth2022_vol_comp,width = 8, height = 5, units = "in", dpi = 300) 

###############################################################
# Edit by Sukhyun Joo ############

# Additional libraries
library(emmeans)
library(nlme)
library(lme4)

unique(growth$event)

# Combine 16-20 and 8-9D
growth$event2 <- growth$event
growth$event2[growth$event == "LC-102 16-20" |
                growth$event == "LC-102 8-9D"] <- "escape"

# Fit model to test a difference between transgenic and control (construct2)
# Check block effect
# Height
mod0 <- lm(H497 ~ H49 + construct2 + block, data = growth)
summary(mod0)

# Make a block as a random effect
mod1 <- lme(H497 ~ H49 + construct2, random = ~1|block, data = growth)
summary(mod1)

# Residual and qq plots
plot(fitted(mod1), residuals(mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod1)); qqline(residuals(mod1))

# Diameter
mod00 <- lm(D497 ~ D49 + construct2 + block, data = growth)
summary(mod00)

# Make a block as a random effect
mod2 <- lme(D497 ~ D49 + construct2, random = ~1|block, data = growth)
summary(mod2)

# Residual and qq plots
plot(fitted(mod2), residuals(mod2), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod2)); qqline(residuals(mod2))

#Volume

mod000 <- lm(log(V497) ~ log(V49) + construct2 + block, data = growth)
summary(mod000)
#Block is significant
exp(0.30454)
mean(growth$V497_cm)
1882.064*0.35

# Make a block as a random effect
mod3 <- lme(log(V497) ~ log(V49) + construct2, random = ~1|block, data = growth)
summary(mod3)

# Residual and qq plots
plot(fitted(mod3), residuals(mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod3)); qqline(residuals(mod3))

##w/o log transformation #####

mod4 <- lme(V497 ~ V49 + construct2, random = ~1|block, data = growth)
summary(mod4)

# Residual and qq plots
plot(fitted(mod4), residuals(mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod4)); qqline(residuals(mod4))



#### Effect of Construct (controls not pooled)####################################
ht_mod3 <- lme(H497 ~ H49 + construct, random = ~1|block, data = growth)
summary(ht_mod3)

# Residual and qq plots
plot(fitted(ht_mod3), residuals(ht_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod3)); qqline(residuals(ht_mod3))

emmeans(ht_mod3, specs = pairwise ~construct)

#effect size

predict_set_h <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","H49","H497"))
predict_set_h$predicted_H497 <- predict(ht_mod3, new_data = predict_set_h, interval = "confidence")

#average modeled height
modelled_control_mean_ht <- mean(predict_set_h$predicted_H497)
modelled_transgenic_mean_ht <- mean(predict_set_h$predicted_H497) + 169.0
modelled_escape_mean_ht <- mean(predict_set_h$predicted_H497) + 92.1

construct_height_effect <- as.data.frame(c(modelled_control_mean_ht, modelled_escape_mean_ht, modelled_transgenic_mean_ht))

row.names(construct_height_effect) <- c("control","escape","transgenic")
colnames(construct_height_effect) <- c("height")


#for interest
ggplot(predict_set_h,aes(x= H497, y = predicted_H497))+
  geom_point()



#diameter construct effect
d_mod3 <- lme(D497 ~ D49 + construct, random = ~1|block, data = growth)
summary(d_mod3)

# Residual and qq plots
plot(fitted(d_mod3), residuals(d_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod3)); qqline(residuals(d_mod3))

emmeans(d_mod3, specs = pairwise ~construct)

predict_set_d <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","D49","D497"))
predict_set_d$predicted_D497 <- predict(d_mod3, new_data = predict_set_d, interval = "confidence")

#average modeled diameter
modelled_control_mean_d <- mean(predict_set_d$predicted_D497)
modelled_transgenic_mean_d <- mean(predict_set_d$predicted_D497) + 1.198
modelled_escape_mean_d <- mean(predict_set_d$predicted_D497) + 1.324

construct_diameter_effect <- as.data.frame(c(modelled_control_mean_d, modelled_escape_mean_d, modelled_transgenic_mean_d))

row.names(construct_diameter_effect) <- c("control","escape","transgenic")
colnames(construct_diameter_effect) <- c("diameter")



###Volume construct effect

v_mod3 <- lme(log(V497) ~ (log(V49)) + construct, random = ~1|block, data = growth)
summary(v_mod3)

# Residual and qq plots
plot(fitted(v_mod3), residuals(v_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod3)); qqline(residuals(v_mod3))

emmeans(v_mod3, specs = pairwise ~construct)
#exp(0.1082)
#exp(0.0820)
#exp(0.0692)
#exp(0.0622)
#exp(0.0390)
#exp(0.0970)

predict_set_v <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V497"))
predict_set_v$predicted_V497 <- exp(predict(v_mod3, new_data = predict_set_v, interval = "confidence"))

#check
ggplot(predict_set_v,aes(x= log(V497), y = predicted_V497))+
  geom_point()
mean(predict_set_v$V497)

#average modeled volume

modelled_control_mean_v <- (mean(predict_set_v$predicted_V497))
modelled_transgenic_mean_v <- ((mean(predict_set_v$predicted_V497))*(exp(0.1082)))
modelled_escape_mean_v <- ((mean(predict_set_v$predicted_V497))*(exp(0.0390)))

construct_volume_effect <- as.data.frame(c(modelled_control_mean_v, modelled_escape_mean_v, modelled_transgenic_mean_v))

row.names(construct_volume_effect) <- c("control","escape","transgenic")
colnames(construct_volume_effect) <- c("volume_index")


######### Effect of event (escapes pooled) ######
ht_mod4 <- lme(H497 ~ H49 + event2, random = ~1|block, data = growth)
summary(ht_mod4)

# Residual and qq plots
plot(fitted(ht_mod4), residuals(ht_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod4)); qqline(residuals(ht_mod4))

emmeans(ht_mod4, specs = pairwise ~event2)

#Effect size
predict_set_h_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","H49","H497"))
predict_set_h_II$predicted_H497 <- predict(ht_mod4, new_data = predict_set_h_II, interval = "confidence")

#average modeled height
modelled_CT3_mean_ht <- mean(predict_set_h_II$predicted_H497)
modelled_pooled_escape_mean_ht <- mean(predict_set_h_II$predicted_H497) + 91.466
modelled_5A_mean_ht <- mean(predict_set_h_II$predicted_H497) + 413.380
modelled_5C_mean_ht <- mean(predict_set_h_II$predicted_H497) + 221.024
modelled_5_mean_ht <- mean(predict_set_h_II$predicted_H497) + 261.446
modelled_2H_mean_ht <- mean(predict_set_h_II$predicted_H497) - 64.454
modelled_13_15E_mean_ht <- mean(predict_set_h_II$predicted_H497) + 32.196
modelled_13_15B_mean_ht <- mean(predict_set_h_II$predicted_H497) + 43.547

event_height_effect <- as.data.frame(c(modelled_CT3_mean_ht, modelled_pooled_escape_mean_ht, modelled_5A_mean_ht, modelled_5C_mean_ht, modelled_5_mean_ht, modelled_2H_mean_ht, modelled_13_15E_mean_ht, modelled_13_15B_mean_ht))

row.names(event_height_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B")
colnames(event_height_effect) <- c("height")

####diameter event effect#######

d_mod4 <- lme(D497 ~ D49 + event2, random = ~1|block, data = growth)
summary(d_mod4)

# Residual and qq plots
plot(fitted(d_mod4), residuals(d_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod4)); qqline(residuals(d_mod4))

emmeans(d_mod4, specs = pairwise ~event2)

#Effect size
predict_set_d_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","D49","D497"))
predict_set_d_II$predicted_D497 <- predict(d_mod4, new_data = predict_set_d_II, interval = "confidence")

#average modeled height
modelled_CT3_mean_d <- mean(predict_set_d_II$predicted_D497)
modelled_pooled_escape_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.3197
modelled_5A_mean_d <- mean(predict_set_d_II$predicted_D497) + 3.0533
modelled_5C_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.9791
modelled_5_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.5744
modelled_2H_mean_d <- mean(predict_set_d_II$predicted_D497) - 0.3084
modelled_13_15E_mean_d <- mean(predict_set_d_II$predicted_D497) + 0.4132
modelled_13_15B_mean_d <- mean(predict_set_d_II$predicted_D497) + 0.5951

event_diameter_effect <- as.data.frame(c(modelled_CT3_mean_d, modelled_pooled_escape_mean_d, modelled_5A_mean_d, modelled_5C_mean_d, modelled_5_mean_d, modelled_2H_mean_d, modelled_13_15E_mean_d, modelled_13_15B_mean_d))

row.names(event_diameter_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B")
colnames(event_diameter_effect) <- c("diameter")


###volume event effect

v_mod4 <- lme(log(V497) ~ log(V49) + event2, random = ~1|block, data=growth)
summary(v_mod4)

plot(fitted(v_mod4), residuals(v_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4)); qqline(residuals(v_mod4))

emmeans(v_mod4, specs = pairwise ~event2)

#5A -CT3
exp(0.26264)
exp(0.1088)
#5A - escape
exp(0.22411)
exp(0.0951)
#5A - 13-15B
exp(0.32854)
exp(0.1249)
#5A - 13-15E
exp(0.28425)
exp(0.1099)
#5A - 2H
exp(0.33052)
exp(0.1107)
#

##Effect size

#Effect size
predict_set_v_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V497"))
predict_set_v_II$predicted_V497 <- exp(predict(v_mod4, new_data = predict_set_v_II, interval = "confidence"))

#check
ggplot(predict_set_v_II,aes(x= log(V497), y = predicted_V497))+
  geom_point()
mean(predict_set_v_II$V497)

#average modeled volume
modelled_CT3_mean_v <- mean(predict_set_v_II$predicted_V497)
modelled_pooled_escape_mean_v <- mean(predict_set_v_II$predicted_V497)*(exp(0.03853))
modelled_5A_mean_v <- mean(predict_set_v_II$predicted_V497)*(exp(0.26264))
modelled_5C_mean_v <- mean(predict_set_v_II$predicted_V497)*(exp(0.16925))
modelled_5_mean_v <- mean(predict_set_v_II$predicted_V497)*(exp(0.23164))
modelled_2H_mean_v <- mean(predict_set_v_II$predicted_V497)/(exp(0.06788))
modelled_13_15E_mean_v <- mean(predict_set_v_II$predicted_V497)/(exp(0.02161))
modelled_13_15B_mean_v <- mean(predict_set_v_II$predicted_V497)/(exp(0.06590))

event_volume_effect <- as.data.frame(c(modelled_CT3_mean_v, modelled_pooled_escape_mean_v, modelled_5A_mean_v, modelled_5C_mean_v, modelled_5_mean_v, modelled_2H_mean_v, modelled_13_15E_mean_v, modelled_13_15B_mean_v))

row.names(event_volume_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B")
colnames(event_volume_effect) <- c("volume")


###############

write.csv(growth,file = "2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")

###table to summarize mean volume per event
Volume_summary <- data_frame()
event_volume <- as.data.frame(aggregate(V497_cm ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(V497_cm ~ event_short, data = growth, FUN = var)
event_n <- aggregate(V497_cm ~ event_short, data = growth, FUN = length)

event_volume_summary <- inner_join(event_volume, event_var, by = "event_short")
event_volume_summary <- inner_join(event_volume_summary, event_n, by = "event_short")

colnames(event_volume_summary) <- c("Event", "mean_Volume","var","n")

event_volume_summary$se <- sqrt(event_volume_summary$var/event_volume_summary$n)
event_volume_summary$tst <- qt(0.975, event_volume_summary$n)
event_volume_summary$sd <- sqrt(event_volume_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_volume_summary$upper_volume <- event_volume_summary$mean_Volume + (event_volume_summary$tst*event_volume_summary$se)
event_volume_summary$lower_volume <- event_volume_summary$mean_Volume - (event_volume_summary$tst*event_volume_summary$se)

write.csv(event_volume_summary, file = "2022_growth&inventory_analylsis/growth_analysis/event_volume_table.csv")

###Volume w/o transformation

# Fit model to test a difference between transgenic and control (construct2)
# Check block effect
# Volume
mod000 <- lm(V497 ~ V49 + construct2 + block, data = growth)
summary(mod000)
#Block is not significant but has influence

# Make a block as a random effect
v_mod3.1 <- lme(V497 ~ V49 + construct2, random = ~1|block, data = growth)
summary(v_mod3.1)

# Residual and qq plots
plot(fitted(v_mod3.1), residuals(v_mod3.1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod3.1)); qqline(residuals(v_mod3.1))

plot(v_mod3.1, which=3)



v_mod4.1 <- lme(V497 ~ V49 + event2, random = ~1|block, data=growth)
summary(v_mod4.1)

plot(fitted(v_mod4.1), residuals(v_mod4.1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4.1)); qqline(residuals(v_mod4.1))
#definitely don't look as good as height and diameter.
#Q-Q plot shows there is a prominent right tail.

emmeans(v_mod4.1, specs = pairwise ~event2)
