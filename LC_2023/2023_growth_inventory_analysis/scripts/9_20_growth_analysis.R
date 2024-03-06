#LC growth analysis update - Sep 21 2023
##This script was used to perform analysis of the Living Carbon trial growth data up to this point in Sep 2023



library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)


###############################################################
# Read in data ############

# Data has been cleaned in 2023_field_season_data_cleaning.R
# Negative growth and missing data was imputted for diameter and height. No outliers were removed
# CT1 was excluded as were trees that were damaged and replaced or died.

#read in data
growth <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
#check data types
growth_check <- growth %>% group_by(event_short)
group_data(growth_check)
str(growth)
growth$D385 <- as.numeric(growth$D385)


###############################################################

#Volume summary stats
sd_vol = sd(growth$V801)
se_vol = (sd(growth$V801))/(sqrt(428))
av_vol = mean(growth$V801) 
CV_vol = (sd_vol/av_vol)*100


###############################################################
# Stats Edit by Sukhyun Joo ############

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
mod0 <- lm(H801 ~ H49 + construct2 + block, data = growth)
summary(mod0)

# Make a block as a random effect
mod1 <- lme(H801 ~ H49 + construct2, random = ~1|block, data = growth)
summary(mod1)
mean(growth$H801)

# Residual and qq plots
plot(fitted(mod1), residuals(mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod1)); qqline(residuals(mod1))

# Diameter
mod00 <- lm(D801 ~ D49 + construct2 + block, data = growth)
summary(mod00)
mean(growth$D801)

# Make a block as a random effect
mod2 <- lme(D801 ~ D49 + construct2, random = ~1|block, data = growth)
summary(mod2)

# Residual and qq plots
plot(fitted(mod2), residuals(mod2), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod2)); qqline(residuals(mod2))

hist(growth$D801)

#Volume

mod000 <- lm(V801 ~ V49 + construct2 + block, data = growth)
summary(mod000)

#Block is significant
# Make a block as a random effect

hist(growth$V801)
#response variable not normal. Skewed to left

#investiate whether log transformation makes senses
#is this skew only seen in some of the data?
#constructs?
ggplot(subset(growth, construct2 == "transgenic"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, construct2 == "control"),aes(x=V801))+
  geom_histogram()
#Appears in both
#Events?
ggplot(subset(growth, event_short == "5A"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, event_short == "CT3"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, event2 == "escape"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, event_short == "4A"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, event_short == "2H"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, event_short == "13-15E"),aes(x=V801))+
  geom_histogram()
#harder to tell with smaller sample sizes

#are standard deviations between groups equal?
growth_event <- growth %>% group_by(event_short) %>% summarize(
  N = n(),
  H_std = sd(H801),
  D_std = sd(D801),
  V_std = sd(V801),
  H = mean(H801),
  D = mean(D801),
  V = mean(V801)
)
#not too bad

#compare max and min values for this response var
max(growth$V801)/min(growth$V801)

#oh wow, yes log transformation necessary

mod3 <- lme(log(V801) ~ V49 + construct2, random = ~1|block, data = growth)
summary(mod3)

# Residual and qq plots
plot(fitted(mod3), residuals(mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod3)); qqline(residuals(mod3))



##w/o log transformation #####

mod4 <- lme(V801 ~ V49 + construct2, random = ~1|block, data = growth)
summary(mod4)

# Residual and qq plots
plot(fitted(mod4), residuals(mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod4)); qqline(residuals(mod4))



#### Effect of Construct (controls not pooled)####################################
ht_mod3 <- lme(H801 ~ H49 + construct, random = ~1|block, data = growth)
summary(ht_mod3)

# Residual and qq plots
plot(fitted(ht_mod3), residuals(ht_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod3)); qqline(residuals(ht_mod3))

emmeans(ht_mod3, specs = pairwise ~construct)

#effect size

predict_set_h <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","H49","H801"))
predict_set_h$predicted_H801 <- predict(ht_mod3, new_data = predict_set_h, interval = "confidence")

#average modeled height
modelled_control_mean_ht <- mean(predict_set_h$predicted_H801)
modelled_transgenic_mean_ht <- mean(predict_set_h$predicted_H801) + 459
modelled_escape_mean_ht <- mean(predict_set_h$predicted_H801) + 338

construct_height_effect <- as.data.frame(c(modelled_control_mean_ht, modelled_escape_mean_ht, modelled_transgenic_mean_ht))

row.names(construct_height_effect) <- c("control","escape","transgenic")
colnames(construct_height_effect) <- c("height")


#diameter construct effect
max(growth$D801)
min(growth$D801)

d_mod3 <- lme(D801 ~ D49 + construct, random = ~1|block, data = growth)
summary(d_mod3)

# Residual and qq plots
plot(fitted(d_mod3), residuals(d_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod3)); qqline(residuals(d_mod3))

emmeans(d_mod3, specs = pairwise ~construct)

predict_set_d <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","D49","D497"))
predict_set_d$predicted_D801 <- predict(d_mod3, new_data = predict_set_d, interval = "confidence")

#average modeled diameter
modelled_control_mean_d <- mean(predict_set_d$predicted_D801)
modelled_transgenic_mean_d <- mean(predict_set_d$predicted_D801) + 2.78
modelled_escape_mean_d <- mean(predict_set_d$predicted_D801) + 4.13

construct_diameter_effect <- as.data.frame(c(modelled_control_mean_d, modelled_escape_mean_d, modelled_transgenic_mean_d))

row.names(construct_diameter_effect) <- c("control","escape","transgenic")
colnames(construct_diameter_effect) <- c("diameter")

###Volume construct effect

v_mod3 <- lme(log(V801) ~ V49 + construct, random = ~1|block, data = growth)
summary(v_mod3)

# Residual and qq plots
plot(fitted(v_mod3), residuals(v_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod3)); qqline(residuals(v_mod3))

v_mod3_emmeans <- emmeans(v_mod3, specs = pairwise ~construct)
v_mod3_emmeans
summary(v_mod3)


#Determine what the effect size difference we saw was...

construct_volume_effect <- v_mod3_emmeans$emmeans
construct_volume_effect <- as.data.frame(v_mod3_emmeans$emmeans)
construct_volume_effect$modelled_vol <- exp(construct_volume_effect$emmean)
construct_volume_effect$SE <- exp(construct_volume_effect$SE)
construct_volume_effect$lower.CL <- exp(construct_volume_effect$lower.CL)
construct_volume_effect$upper.CL <- exp(construct_volume_effect$upper.CL)


#interpreting the power of the construct comparisons

#power in comparing between transgenic and escape?
power_LC102_escape <- subset(growth, construct == "LC-102"|construct =="Escape") %>% group_by(construct) %>% summarize(
  n = n(),
  var = var(V801),
  mean = mean(V801),
  sd = sd(V801),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

LC102_escape_diff_mean <- (power_LC102_escape[2,4] - power_LC102_escape[1,4])
LC102_escape_diff_lower <- (power_LC102_escape[2,8] - power_LC102_escape[1,8])
LC102_escape_diff_upper <- (power_LC102_escape[2,9] - power_LC102_escape[1,9])
pooled_var_LC_102_escape <- ((power_LC102_escape[1,2]*power_LC102_escape[1,3]) + (power_LC102_escape[2,2]*power_LC102_escape[2,3]))/(power_LC102_escape[1,2]+power_LC102_escape[2,2]-2)


sig_effect_size_LC_102_Escape <- sqrt((pooled_var_LC_102_escape[1,1])*2/62)*(qt(0.975,61))/(mean(c(as.numeric(power_LC102_escape[1,4]),as.numeric(power_LC102_escape[2,4]))))

#we had the power to detect an effect size of 28.8%

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))
sample_calc <- function(df,df2,df3) {
  av = mean(c(as.numeric(df[1,4]),as.numeric(df[2,4])))
  effect_5 = (av*(df2[1,1]/100))
  effect_10 = (av*(df2[2,1]/100))
  effect_15 = (av*(df2[3,1]/100))
  effect_20 = (av*(df2[4,1]/100))
  effect_25 = (av*(df2[5,1]/100))
  effect_30 = av*(df2[6,1]/100)
  effect_35 = av*(df2[7,1]/100)
  effect_40 = av*(df2[8,1]/100)
  effect_45 = av*(df2[9,1]/100)
  effect_50 = av*(df2[10,1]/100)
  n_5 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_5)^2)
  n_10 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_10)^2)
  n_15 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_15)^2)
  n_20 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_20)^2)
  n_25 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_25)^2)
  n_30 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_30)^2)
  n_35 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_35)^2)
  n_40 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_40)^2)
  n_45 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_45)^2)
  n_50 = ((qt(0.975,61)*sqrt(2*(df3[1,1]))/effect_50)^2)
  effect_size = c(effect_5,effect_10, effect_15, effect_20, effect_25, effect_30, effect_35,effect_40, effect_45, effect_50)
  sample_size = c(n_5,n_10,n_15,n_20,n_25,n_30,n_35,n_40,n_45,n_50)
  effect_table_df =data.frame(df2,effect_percent,sample_size)
  colnames(effect_table_df) = c("percent","effect","n")
  return(effect_table_df)
}
LC102_escape_power_table <- sample_calc(power_LC102_escape, effect_percent,pooled_var_LC_102_escape)

ggplot(subset(LC102_escape_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1/µ1) (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=62),color = "blue",linetype = 2)+
  geom_line(aes(x=10.7),color = "red",linetype =2)


######### Effect of event (escapes pooled) ######
ht_mod4 <- lme(H801 ~ H49 + event2, random = ~1|block, data = growth)
summary(ht_mod4)

# Residual and qq plots
plot(fitted(ht_mod4), residuals(ht_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod4)); qqline(residuals(ht_mod4))

ht_mod4_emmeans <- emmeans(ht_mod4, specs = pairwise ~event2)
ht_mod4_emmeans
# effect size
event_height_effect <- as.data.frame(ht_mod4_emmeans$emmeans)
event_height_effect$modelled_ht <- event_height_effect$emmean


####diameter event effect#######

d_mod4 <- lme(D801 ~ D49 + event2, random = ~1|block, data = growth)
summary(d_mod4)

hist(growth$D801)
# Residual and qq plots
plot(fitted(d_mod4), residuals(d_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod4)); qqline(residuals(d_mod4))

d_mod4_emmeans <- emmeans(d_mod4, specs = pairwise ~event2)
d_mod4_emmeans

# effect size
event_diameter_effect <- as.data.frame(d_mod4_emmeans$emmeans)
event_diameter_effect$modelled_ht <- event_diameter_effect$emmean


###volume event effect

v_mod4 <- lme(log(V801) ~ V49 + event2, random = ~1|block, data=growth)
summary(v_mod4)

plot(fitted(v_mod4), residuals(v_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4)); qqline(residuals(v_mod4))

v_mod4_emmeans <- emmeans(v_mod4, specs = pairwise ~event2)
v_mod4_emmeans

#Effect size
predict_set_v_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V801"))
predict_set_v_II$predicted_V801 <- exp(predict(v_mod4, new_data = predict_set_v_II, interval = "confidence"))

# Combine 16-20 and 8-9D
predict_set_v_II$event2 <- predict_set_v_II$event
predict_set_v_II$event2[predict_set_v_II$event == "LC-102 16-20" |
                          predict_set_v_II$event == "LC-102 8-9D"] <- "escape"

#average modeled volume
event_volume_effect <- predict_set_v_II %>% group_by(event2) %>% summarize(
  predicted_mean = mean(predicted_V801),
  mean = mean(V801)
)

event_volume_effect <- as.data.frame(v_mod4_emmeans$emmeans)
event_volume_effect$modelled_vol <- exp(event_volume_effect$emmean)
event_volume_effect$SE <- exp(event_volume_effect$SE)
event_volume_effect$lower.CL <- exp(event_volume_effect$lower.CL)
event_volume_effect$upper.CL <- exp(event_volume_effect$upper.CL)

#With our sampling sizes, what difference in volume were we able to detect as statistically different from 0?

#power in comparing between 5A and escape
######what is the power we had when we used our model to account for variation in starting volume and block

power_5A_escape <- subset(growth, event2 == "LC-102 5A"|event2 =="escape") %>% group_by(event2) %>% summarize(
  n = n(),
  var = var(V801),
  mean = mean(V801),
  sd = sd(V801),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

escape_5A_diff_mean <- (power_5A_escape[1,4] - power_5A_escape[2,4])
escape_5A_diff_lower <- (power_5A_escape[1,8] - power_5A_escape[2,8])
escape_5A_diff_upper <- (power_5A_escape[1,9] - power_5A_escape[2,9])
pooled_var_5A_escape <- ((power_5A_escape[1,2]*power_5A_escape[1,3]) + (power_5A_escape[2,2]*power_5A_escape[2,3]))/(power_5A_escape[1,2]+power_5A_escape[2,2]-2)


sig_effect_size_5A_Escape <- sqrt(((pooled_var_5A_escape[1,1])*2)/34)*(qt(0.975,33))/(mean(c(as.numeric(power_5A_escape[1,4]),as.numeric(power_5A_escape[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

escape_5A_power_table2 <- sample_calc(power_5A_escape, effect_percent,pooled_var_5A_escape)

ggplot(subset(escape_5A_power_table2, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=34), color="blue",linetype = 2)+
  geom_line(aes(x=14), color = "red", linetype = 2)


############# 5A and CT3 ##################

power_5A_CT3 <- subset(growth, event2 == "LC-102 5A"|event2 =="CT717 3") %>% group_by(event2) %>% summarize(
  n = n(),
  var = var(V801),
  mean = mean(V801),
  sd = sd(V801),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

CT3_5A_diff_mean <- (power_5A_CT3[2,4] - power_5A_CT3[1,4])
CT3_5A_diff_lower <- (power_5A_CT3[2,8] - power_5A_CT3[1,8])
CT3_5A_diff_upper <- (power_5A_CT3[2,9] - power_5A_CT3[1,9])
pooled_var_5A_CT3 <- ((power_5A_CT3[1,2]*power_5A_CT3[1,3]) + (power_5A_CT3[2,2]*power_5A_CT3[2,3]))/(power_5A_CT3[1,2]+power_5A_CT3[2,2]-2)


sig_effect_size_5A_CT3 <- sqrt((pooled_var_5A_CT3[1,1])*2/33)*(qt(0.975,32))/(mean(c(as.numeric(power_5A_CT3[1,4]),as.numeric(power_5A_CT3[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

CT3_5A_power_table <- sample_calc(power_5A_CT3, effect_percent,pooled_var_5A_CT3)

ggplot(subset(CT3_5A_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  
  theme_bw()+
  geom_line(aes(y=34), color="blue",linetype = 2)+
  geom_line(aes(x=14), color = "red", linetype = 2)


##############################
#####comparing event tiers or classes or groupings or whatever
growth <- growth %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D"~"control",
  event_short == "CT3"~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))

growth$tier <- as.factor(growth$tier)

########height########
ht_mod5 <- lme(H801 ~ H49 + tier, random = ~1|block, data = growth)
summary(ht_mod5)

# Residual and qq plots
plot(fitted(ht_mod5), residuals(ht_mod5), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod5)); qqline(residuals(ht_mod5))

ht_mod5_emmeans <- emmeans(ht_mod5, specs = pairwise ~tier)
ht_mod5_emmeans
# effect size
tier_height_effect <- as.data.frame(ht_mod5_emmeans$emmeans)
tier_height_effect$modelled_ht <- tier_height_effect$emmean

###diameter
d_mod5 <- lme(D801 ~ D49 + tier, random = ~1|block, data = growth)
summary(d_mod5)

# Residual and qq plots
plot(fitted(d_mod5), residuals(d_mod5), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod5)); qqline(residuals(d_mod5))

d_mod5_emmeans <- emmeans(d_mod5, specs = pairwise ~tier)
d_mod5_emmeans
# effect size
tier_diameter_effect <- as.data.frame(d_mod5_emmeans$emmeans)
tier_diameter_effect$modelled_diameter <- tier_diameter_effect$emmean


########################volume######

v_mod5 <- lme(log(V801) ~ V49 + tier, random = ~1|block, data = growth)
summary(v_mod5)

# Residual and qq plots
plot(fitted(v_mod5), residuals(v_mod5), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod5)); qqline(residuals(v_mod5))

v_mod5_emmeans <- emmeans(v_mod5, specs = pairwise ~tier)
v_mod5_emmeans
# effect size
tier_volume_effect <- as.data.frame(v_mod5_emmeans$emmeans)
tier_volume_effect$modelled_volume <- exp(tier_volume_effect$emmean)

####what effect size would have been needed to observe sig results when comparing top events to escapes

power_top_escape <- subset(growth, tier == "top"|tier =="control") %>% group_by(tier) %>% summarize(
  n = n(),
  var = var(V801),
  mean = mean(V801),
  sd = sd(V801),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

top_escape_diff_mean <- (power_top_escape[2,4] - power_top_escape[1,4])
top_escape_diff_lower <- (power_top_escape[2,8] - power_top_escape[1,8])
top_escape_diff_upper <- (power_top_escape[2,9] - power_top_escape[1,9])
pooled_var_top_escape <- ((power_top_escape[1,2]*power_top_escape[1,3]) + (power_top_escape[2,2]*power_top_escape[2,3]))/(power_top_escape[1,2]+power_top_escape[2,2]-2)


sig_effect_size_top_escape <- sqrt((pooled_var_top_escape[1,1])*2/62)*(qt(0.975,61))/(mean(c(as.numeric(power_top_escape[1,4]),as.numeric(power_top_escape[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

top_escape_power_table <- sample_calc(power_top_escape, effect_percent,pooled_var_top_escape)

top_escape_power_plot <- ggplot(subset(top_escape_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=62), color="blue",linetype = 2)+
  geom_line(aes(x=9.4), color = "red", linetype = 2)

top_escape_power_plot
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/top_escape_power_plot.png",plot=top_escape_power_plot, dpi = 300)

#########################################
#partition variance between tier, and environment for ht, diameter and volume
#Height
aov0 <- aov(H801 ~ H49 + tier + block, data = growth)
anova(aov0)
var(growth$H801)
#Diameter
aov00 <- aov(D801 ~ D49 + tier + block, data = growth)
anova(aov00)

#Volume
aov000 <- aov(V801 ~ V49 + tier  + block, data = growth)
anova(aov000)

##########################################

#yearly increments
growth$yr1_vol <- growth$V299- growth$V49
growth$yr2_vol <- growth$V497 - growth$V299
growth$yr3_vol <- growth$V801 - growth$V497

volume_yearly_summary <- growth %>% group_by(event_short) %>% dplyr::summarise(
  n = n(),
  yr1_vol_sd = sd(yr1_vol),
  yr2_vol_sd = sd(yr2_vol),
  yr3_vol_sd = sd(yr3_vol),
  yr1_vol = mean(yr1_vol),
  yr2_vol = mean(yr2_vol),
  yr3_vol = mean(yr3_vol),
  yr1_vol_se = yr1_vol_sd/(sqrt(n)),
  yr2_vol_se = yr2_vol_sd/(sqrt(n)),
  yr3_vol_se = yr3_vol_sd/(sqrt(n))
)

volume_yearly_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

myColors <- c('gray2',
              'cyan4',
              'red4')


names(myColors) <- levels(volume_yearly_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)

## growth year 1
growth_year1 <- ggplot(volume_yearly_summary, aes(x=reorder(event_short,yr1_vol),y=yr1_vol, color = construct))+
  geom_point()+
  geom_errorbar(aes(ymin=yr1_vol-(2*yr1_vol_se), ymax=yr1_vol+(2*yr1_vol_se)))+
  xlab("event")+
  ylab("Year one volume index (cm3)")

growth_year1 + colorScale
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/growth_yr1.png",plot=growth_year1+colorScale, dpi=300)

##growth year 2
growth_year2 <- ggplot(volume_yearly_summary, aes(x=reorder(event_short,yr2_vol),y=yr2_vol, color = construct))+
  geom_point()+
  xlab("event")+
  ylab("Year one volume index (cm3)")+
  geom_errorbar(aes(ymin=yr2_vol-(2*yr2_vol_se), ymax=yr2_vol+(2*yr2_vol_se)))

growth_year2 + colorScale
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/growth_yr2.png",plot=growth_year2+colorScale, dpi=300)

#growth year 3
growth_year3 <- ggplot(volume_yearly_summary, aes(x=reorder(event_short,yr3_vol),y=yr3_vol, color = construct))+
  geom_point()+
  geom_errorbar(aes(ymin=yr3_vol-(2*yr3_vol_se), ymax=yr3_vol+(2*yr3_vol_se)))+
  xlab("event")+
  ylab("Year one volume index (cm3)")

growth_year3 + colorScale
ggsave(filename = "LC_2023/2023_growth_inventory_analysis/growth_yr3.png",plot=growth_year3+colorScale, dpi=300)

######
#Height to diameter ratio
growth$HD801 <- growth$H801/growth$D801
growth$HD49 <- growth$H49/growth$V49

HD_mod5 <- lme(HD801 ~ HD49 + tier, random = ~1|block, data = growth)
summary(HD_mod5)

# Residual and qq plots
plot(fitted(HD_mod5), residuals(HD_mod5), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(HD_mod5)); qqline(residuals(HD_mod5))

HD_mod5_emmeans <- emmeans(HD_mod5, specs = pairwise ~tier)
HD_mod5_emmeans

###event comparison


HD_mod4 <- lme(HD801 ~ HD49 + event2, random = ~1|block, data = growth)
summary(HD_mod4)

hist(growth$HD801)
# Residual and qq plots
plot(fitted(HD_mod4), residuals(HD_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(HD_mod4)); qqline(residuals(HD_mod4))

HD_mod4_emmeans <- emmeans(HD_mod4, specs = pairwise ~event2)
HD_mod4_emmeans

# effect size
event_HD_effect <- as.data.frame(HD_mod4_emmeans$emmeans)
event_HD_effect$modelled_HD <- event_HD_effect$emmean


