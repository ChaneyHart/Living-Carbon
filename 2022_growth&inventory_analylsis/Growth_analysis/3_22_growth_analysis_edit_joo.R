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



ht_full <- ggplot(Height_long, aes(x = Days, y = meters))+
  geom_point(aes(color = ID), show.legend = FALSE)+
  xlab("Days since planting")+
  ylab("Height (m)")


ht_full
ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline.png", plot = ht_full, width = 6, height = 3, units = "in", dpi = 300)

d_full <- ggplot(Diam_long, aes(x = Days, y = Diameter))+
  geom_point(aes(color = ID), show.legend = FALSE)+
  xlab("days since planting")+
  ylab("Diameter (mm)")


d_full
ggsave("2022_growth&inventory_analylsis/growth_graph/d_full.png", plot = d_full, width = 6, height = 3, units = "in", dpi = 300)


vol_full <- ggplot(Volume_long, aes(x = Days, y = log(Volume)))+
  geom_point(aes(color = ID), show.legend = FALSE)+
  xlab("days since planting")+
  ylab("Volume log(cubic mm)")



vol_full
ggsave("2022_growth&inventory_analylsis/growth_graph/vol_full.png", plot = vol_full, width = 6, height = 3, units = "in", dpi = 300)


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

# Make a block as a random effect
mod3 <- lme(log(V497) ~ log(V49) + construct2, random = ~1|block, data = growth)
summary(mod3)

# Residual and qq plots
plot(fitted(mod3), residuals(mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod3)); qqline(residuals(mod3))
#

# Effect of Construct (controls not pooled)
ht_mod3 <- lme(H497 ~ H49 + construct, random = ~1|block, data = growth)
summary(ht_mod3)

emmeans(ht_mod3, specs = pairwise ~construct)

# Residual and qq plots
plot(fitted(ht_mod3), residuals(ht_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod3)); qqline(residuals(ht_mod3))

d_mod3 <- lme(D497 ~ D49 + construct, random = ~1|block, data = growth)
summary(d_mod3)

emmeans(d_mod3, specs = pairwise ~construct)


# Residual and qq plots
plot(fitted(d_mod3), residuals(d_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod3)); qqline(residuals(d_mod3))


v_mod3 <- lme((V497) ~ (V49) + construct, random = ~1|block, data = growth)
summary(v_mod3)

emmeans(v_mod3, specs = pairwise ~construct)
#exp(0.1082)
#exp(0.0820)
#exp(0.0692)
#exp(0.0622)
#exp(0.0390)
#exp(0.0970)


# Effect of event (escapes pooled)
ht_mod4 <- lme(H497 ~ H49 + event2, random = ~1|block, data = growth)
summary(ht_mod4)

# Residual and qq plots
plot(fitted(ht_mod4), residuals(ht_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod4)); qqline(residuals(ht_mod4))

emmeans(ht_mod4, specs = pairwise ~event2)

d_mod4 <- lme(D497 ~ D49 + event2, random = ~1|block, data = growth)
summary(d_mod4)

# Residual and qq plots
plot(fitted(d_mod4), residuals(d_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod4)); qqline(residuals(d_mod4))

emmeans(d_mod4, specs = pairwise ~event2)

v_mod4 <- lme(log(V497) ~ log(V49) + event2, random = ~1|block, data=growth)
summary(v_mod4)

plot(fitted(v_mod4), residuals(v_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4)); qqline(residuals(v_mod4))
#definitely don't look as good as height and diameter.
#Q-Q plot shows there is a prominent right tail.

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

write.csv(growth,file = "2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")

# Fit model to test a difference between transgenic and control (construct2)
# Check block effect
# Volume
mod000 <- lm(log(V497) ~ log(V49) + construct2 + block, data = growth)
summary(mod000)
#Block is significant

# Make a block as a random effect
v_mod3 <- lme(V497 ~ V49 + construct2, random = ~1|block, data = growth)
summary(v_mod3)

# Residual and qq plots
plot(fitted(mod3), residuals(mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod3)); qqline(residuals(mod3))
#looks ok

