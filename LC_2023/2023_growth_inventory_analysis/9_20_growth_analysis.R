#LC growth analysis update - Sep 21 2023
##This script was used to perform analysis of the Living Carbon trial growth data up to this point in Sep 2023



library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)
library(multcomp)



#bring in growth data w/ info about each tree
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

#calculate volume


growth$V49 = ((3.14159*(((growth$D49/10)/2)^2))*(growth$H49/10))
growth$V144 = ((3.14159*(((growth$D144/10)/2)^2))*(growth$H144/10))
growth$V299 = ((3.14159*(((growth$D299/10)/2)^2))*(growth$H299/10))
growth$V335 = ((3.14159*(((growth$D335/10)/2)^2))*(growth$H335/10))
growth$V357 = ((3.14159*(((growth$D357/10)/2)^2))*(growth$H357/10))
growth$V385 = ((3.14159*(((growth$D385/10)/2)^2))*(growth$H385/10))
growth$V421 = ((3.14159*(((growth$D421/10)/2)^2))*(growth$H421/10))
growth$V497 = ((3.14159*(((growth$D497/10)/2)^2))*(growth$H497/10))
growth$V664 = ((3.14159*(((growth$D664/10)/2)^2))*(growth$H664/10))
growth$V801 = ((3.14159*(((growth$D801/10)/2)^2))*(growth$H801/10))



#Convert dataframe to long format to look at changes over time for each tree
Height_long <- pivot_longer(growth, cols = c(H49,H144,H299,H335,H357,H385,H421,H497,H664,H801), names_to = "Days", values_to = "Height")
Height_long <- Height_long %>% mutate(Days = as.numeric(substr(Height_long$Days,2,4)))
Height_long$meters <- Height_long$Height/1000

Diam_long <- pivot_longer(growth, cols = c(D49,D144,D299,D335,D357,D385,D421,D497,D664,D801), names_to = "Days", values_to = "Diameter")
Diam_long <- Diam_long %>% mutate(Days = as.numeric(substr(Diam_long$Days,2,4)))

Volume_long <- pivot_longer(growth, cols = c(V49,V144,V299,V335,V357,V385,V421,V497,V664,V801), names_to = "Days", values_to = "Volume")
Volume_long <- Volume_long %>% mutate(Days = as.numeric(substr(Volume_long$Days,2,4)))


#Create a custom color scale
  library(RColorBrewer)
myColors <- c('cyan4',
              'red4',
              'gray')


names(myColors) <- levels(growth$construct)
colorScale <- scale_fill_manual(name = "construct",values = myColors)

ht_full <- ggplot(Height_long, aes(x = Days, y = meters))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Height (m)")


ht_full
ggsave("LC_2023/2023_growth_inventory_analysis/ht_timeline_2023.png", plot = ht_full, width = 6, height = 3, units = "in", dpi = 300)

d_full <- ggplot(Diam_long, aes(x = Days, y = Diameter))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("days since planting")+
  ylab("Diameter (mm)")


d_full
ggsave("LC_2023/2023_growth_inventory_analysis/diam_timeline_2023.png", plot = d_full, width = 6, height = 3, units = "in", dpi = 300)


vol_full <- ggplot(Volume_long, aes(x = Days, y = (Volume)))+
  geom_point(aes(color = event_short), show.legend = FALSE, size = 0.75)+
  xlab("days since planting")+
  ylab("Volume (cubic cm)")



vol_full
ggsave("LC_2023/2023_growth_inventory_analysis/vol_timeline_2023.png", plot = vol_full, width = 6, height = 3, units = "in", dpi = 300)

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
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Height (m)")
#large block
ht_full_event_lb <- ggplot(subset(Height_long_II_summary, block == "large"), aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("height (m)")
##small block
ht_full_event_sb <- ggplot(subset(Height_long_II_summary, block == "small"), aes(x = Days, y = height))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("height (m)")



ht_full_event
ht_full_event_lb
ht_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis//ht_timeline_event_2023.png", plot = ht_full_event, width = 6, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_lb.png", plot = ht_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/ht_timeline_event_sb.png", plot = ht_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)



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
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

#large block
Diam_full_event_lb <- ggplot(subset(diam_long_II_summary, block == "large"), aes(x = Days, y = diam))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")
##small block
Diam_full_event_sb <- ggplot(subset(diam_long_II_summary, block == "small"), aes(x = Days, y = diam))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("diameter (mm)")

Diam_full_event
Diam_full_event_lb
Diam_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis//Diam_timeline_event_2023.png", plot = Diam_full_event, width = 6, height = 3, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_lb.png", plot = Diam_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/Diam_timeline_event_sb.png", plot = Diam_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)


#####volume by event ###########

volume_long_II <- Volume_long %>% group_by(event_short,Days, block)
volume_long_III <- Volume_long %>% group_by(event_short,Days)

volume_long_II_summary <- volume_long_II %>% dplyr::summarise(
  vol = mean(Volume),
  n = n(),
  vol_sd = sd(Volume),
  vol_se = vol_sd/(sqrt(n))
)

volume_long_III_summary <- volume_long_III %>% dplyr::summarise(
  vol = mean(Volume),
  n = n(),
  vol_sd = sd(Volume),
  vol_se = vol_sd/(sqrt(n))
)


vol_full_event <- ggplot(volume_long_III_summary, aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")+
  theme_bw()

vol_full_event
#large block
vol_full_event_lb <- ggplot(subset(volume_long_II_summary, block == "large"), aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")

#small block

vol_full_event_sb <- ggplot(subset(volume_long_II_summary, block == "small"), aes(x = Days, y = vol))+
  geom_point(aes(color = event_short), show.legend = TRUE, size = 0.75)+
  xlab("Days since planting")+
  ylab("Volume index (cubic cm)")


vol_full_event
vol_full_event_lb
vol_full_event_sb
ggsave("LC_2023/2023_growth_inventory_analysis/vol_timeline_event_2023.png", plot = vol_full_event, width = 6, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_lb.png", plot = vol_full_event_lb, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("2022_growth&inventory_analylsis/growth_graph/vol_timeline_event_sb.png", plot = vol_full_event_sb, width = 5, height = 4, units = "in", dpi = 300)



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
mean(growth$V497)
V497_se <- sqrt(var(growth$V497)/(length(growth$V497)))
V299_se <- sqrt(var(growth$V299)/(length(growth$V299)))



v_percent_22_upper <- (((mean(growth$V497_cm)+(2*V497_se)) - (mean(growth$V299_cm)+(2*V299_se)))/(mean(growth$V299_cm)+(2*V299_se)))*100
v_percent_22_lower <- (((mean(growth$V497_cm)-(2*V497_se)) - (mean(growth$V299_cm)-(2*V299_se)))/(mean(growth$V299_cm)-(2*V299_se)))*100
v_percent_22 <- (((mean(growth$V497)) - (mean(growth$V299)))/(mean(growth$V299)))*100



(median(growth$V497_cm) - median(growth$V299_cm))/(median(growth$V299_cm))

mean(growth$V299)
mean(growth$V497)
##2023

#2023 height growth
mean(growth$V801)
V801_se <- sqrt(var(growth$V801)/(length(growth$V801)))

ht_percent_23_upper <- (((mean(growth$V801)+(2*H801_se)) - (mean(growth$V497)+(2*V497_se)))/(mean(growth$V497)+(2*V497_se)))*100
ht_percent_23_lower <- (((mean(growth$v801)-(2*H801_se)) - (mean(growth$H497)-(2*H497_se)))/(mean(growth$H497)-(2*H497_se)))*100
vol_percent_23 <- (((mean(growth$V801)) - (mean(growth$V497)))/(mean(growth$V497)))*100

###statistical analysis##############



#Are there any events that grew significantly more in any one year?

str(growth)

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
  geom_errorbar(aes(ymin=yr1_vol-(2*yr1_vol_se), ymax=yr1_vol+(2*yr1_vol_se)))

growth_year1 + colorScale

##growth year 2
growth_year2 <- ggplot(volume_yearly_summary, aes(x=reorder(event_short,yr2_vol),y=yr2_vol, color = construct))+
  geom_point()+
  geom_errorbar(aes(ymin=yr2_vol-(2*yr2_vol_se), ymax=yr2_vol+(2*yr2_vol_se)))

growth_year2 + colorScale

#growth year 3
growth_year3 <- ggplot(volume_yearly_summary, aes(x=reorder(event_short,yr3_vol),y=yr3_vol, color = construct))+
  geom_point()+
  geom_errorbar(aes(ymin=yr3_vol-(2*yr3_vol_se), ymax=yr3_vol+(2*yr3_vol_se)))

growth_year3 + colorScale

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

# Residual and qq plots
plot(fitted(mod1), residuals(mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod1)); qqline(residuals(mod1))

# Diameter
mod00 <- lm(D801 ~ D49 + construct2 + block, data = growth)
summary(mod00)

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
#summary stats for volume 
sd_vol = sd(growth$V801)
se_vol = (sd(growth$V801))/(sqrt(428))
av_vol = mean(growth$V801) 
CV_vol = (sd_vol/av_vol)*100


mod000 <- lm(V801 ~ V49 + construct2 + block, data = growth)
summary(mod000)
#Block is significant


# Make a block as a random effect

hist(growth$V801)
#response variable not normal. Skewed to left

#is this skew only seen in some of the data?
ggplot(subset(growth, construct2 == "transgenic"),aes(x=V801))+
  geom_histogram()
ggplot(subset(growth, construct2 == "control"),aes(x=V801))+
  geom_histogram()


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

emmeans(v_mod3, specs = pairwise ~construct)
summary(v_mod3)

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

ggplot(subset(LC102_escape_power_table, effect > 10), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Effect size (%)")+
  ylab("Sample size needed")+
  theme_bw()
  

#Determine what the effect size difference we saw was...

predict_set_v <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V801"))
predict_set_v$predicted_V801 <- exp(predict(v_mod3, new_data = predict_set_v, interval = "confidence"))


#average modeled volume

modelled_control_mean_v <- (mean(predict_set_v$predicted_V801))
modelled_transgenic_mean_v <- ((mean(predict_set_v$predicted_V801))*(exp(0.234)))
modelled_escape_mean_v <- ((mean(predict_set_v$predicted_V801))*(exp(0.336)))

construct_volume_effect <- as.data.frame(c(modelled_control_mean_v, modelled_escape_mean_v, modelled_transgenic_mean_v))

row.names(construct_volume_effect) <- c("control","escape","transgenic")
colnames(construct_volume_effect) <- c("volume_index")



######### Effect of event (escapes pooled) ######
ht_mod4 <- lme(H801 ~ H49 + event2, random = ~1|block, data = growth)
summary(ht_mod4)

# Residual and qq plots
plot(fitted(ht_mod4), residuals(ht_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(ht_mod4)); qqline(residuals(ht_mod4))

emmeans(ht_mod4, specs = pairwise ~event2)

#Effect size
predict_set_h_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","H49","H801"))
predict_set_h_II$predicted_H801 <- predict(ht_mod4, new_data = predict_set_h_II, interval = "confidence")

#average modeled height
modelled_CT3_mean_ht <- mean(predict_set_h_II$predicted_H801)
modelled_pooled_escape_mean_ht <- mean(predict_set_h_II$predicted_H801) + 335.44
modelled_5A_mean_ht <- mean(predict_set_h_II$predicted_H801) + 859.93
modelled_5C_mean_ht <- mean(predict_set_h_II$predicted_H801) + 565.46
modelled_5_mean_ht <- mean(predict_set_h_II$predicted_H801) + 524.55
modelled_4A_mean_ht <- mean(predict_set_h_II$predicted_H801) + 703.44
modelled_1_mean_ht <- mean(predict_set_h_II$predicted_H801) + 352.9
modelled_1C_mean_ht <- mean(predict_set_h_II$predicted_H801) + 430.96
modelled_4B_mean_ht <- mean(predict_set_h_II$predicted_H801) + 497.24
modelled_7_mean_ht <- mean(predict_set_h_II$predicted_H801) + 488.71
modelled_2H_mean_ht <- mean(predict_set_h_II$predicted_H801) + 80.19
modelled_13_15E_mean_ht <- mean(predict_set_h_II$predicted_H801) + 229.46
modelled_13_15B_mean_ht <- mean(predict_set_h_II$predicted_H801) + 267.59

event_height_effect <- as.data.frame(c(modelled_CT3_mean_ht, modelled_pooled_escape_mean_ht, modelled_5A_mean_ht, modelled_5C_mean_ht, modelled_5_mean_ht, modelled_2H_mean_ht, modelled_13_15E_mean_ht, modelled_13_15B_mean_ht, modelled_1_mean_ht, modelled_1C_mean_ht, modelled_4A_mean_ht, modelled_4B_mean_ht, modelled_7_mean_ht))

row.names(event_height_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B","1","1C","4A","4B","7")
colnames(event_height_effect) <- c("height")

####diameter event effect#######

d_mod4 <- lme(D801 ~ D49 + event2, random = ~1|block, data = growth)
summary(d_mod4)

hist(growth$D801)
# Residual and qq plots
plot(fitted(d_mod4), residuals(d_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod4)); qqline(residuals(d_mod4))

emmeans(d_mod4, specs = pairwise ~event2)

#Effect size
predict_set_d_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","D49","D801"))
predict_set_d_II$predicted_D801 <- predict(d_mod4, new_data = predict_set_d_II, interval = "confidence")

ggplot(predict_set_d_II, aes(x=predicted_D801, y = (D801)))+
  geom_point()

#average modeled height
modelled_CT3_mean_d <- mean(predict_set_d_II$predicted_D801)
modelled_pooled_escape_mean_d <- mean(predict_set_d_II$predicted_D801) + 4.1043
modelled_5A_mean_d <- mean(predict_set_d_II$predicted_D801) + 6.3998
modelled_5C_mean_d <- mean(predict_set_d_II$predicted_D801) + 3.7994
modelled_5_mean_d <- mean(predict_set_d_II$predicted_D801) + 2.5071
modelled_4A_mean_d <- mean(predict_set_d_II$predicted_D801) + 4.0694
modelled_2H_mean_d <- mean(predict_set_d_II$predicted_D801) + 0.5401
modelled_13_15E_mean_d <- mean(predict_set_d_II$predicted_D801) + 1.4315
modelled_13_15B_mean_d <- mean(predict_set_d_II$predicted_D801) + 0.6725 
modelled_1_mean_d <- mean(predict_set_d_II$predicted_D801) + 1.6158
modelled_1C_mean_d <- mean(predict_set_d_II$predicted_D801) + 1.9290
modelled_4B_mean_d <- mean(predict_set_d_II$predicted_D801) + 3.4374
modelled_7_mean_d <- mean(predict_set_d_II$predicted_D801) + 3.2444

event_diameter_effect <- as.data.frame(c(modelled_CT3_mean_d, modelled_pooled_escape_mean_d, modelled_5A_mean_d, modelled_5C_mean_d, modelled_5_mean_d, modelled_4A_mean_d,modelled_2H_mean_d, modelled_13_15E_mean_d, modelled_13_15B_mean_d, modelled_1_mean_d, modelled_1C_mean_d, modelled_4B_mean_d, modelled_7_mean_d))

row.names(event_diameter_effect) <- c("control","escape","5A","5C","5","4A","2H","13_15E","13_15B","1","1C","4B","7")
colnames(event_diameter_effect) <- c("diameter")


###volume event effect

v_mod4 <- lme(log(V801) ~ V49 + event2, random = ~1|block, data=growth)
summary(v_mod4)

plot(fitted(v_mod4), residuals(v_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4)); qqline(residuals(v_mod4))

emmeans(v_mod4, specs = pairwise ~event2)

#Effect size
predict_set_v_II <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V801"))
predict_set_v_II$predicted_V801 <- exp(predict(v_mod4, new_data = predict_set_v_II, interval = "confidence"))


#average modeled volume
modelled_CT3_mean_v <- mean(predict_set_v_II$predicted_V801)
modelled_pooled_escape_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.33449))
modelled_5A_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.47141))
modelled_5C_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.38730))
modelled_5_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.24145))
modelled_2H_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.03951))
modelled_13_15E_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.01260))
modelled_13_15B_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.03764))
modelled_4A_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.41037))
modelled_4B_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.39937))
modelled_1_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.14146))
modelled_1C_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.18226))
modelled_7_mean_v <- mean(predict_set_v_II$predicted_V801)*(exp(0.17645))

event_volume_effect <- as.data.frame(c(modelled_CT3_mean_v, modelled_pooled_escape_mean_v, modelled_5A_mean_v, modelled_5C_mean_v, modelled_5_mean_v, modelled_2H_mean_v, modelled_13_15E_mean_v, modelled_13_15B_mean_v, modelled_4A_mean_v, modelled_4B_mean_v, modelled_1_mean_v, modelled_1C_mean_v, modelled_7_mean_v))

row.names(event_volume_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B","4A","4B","1","1C","7")
colnames(event_volume_effect) <- c("volume")

#With our sampling sizes, what difference in volume were we able to detect as statistically different from 0?

#power in comparing between 5A and escape?

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


sig_effect_size_5A_Escape <- sqrt((pooled_var_5A_escape[1,1])*2/34)*(qt(0.975,33))/(mean(c(as.numeric(power_5A_escape[1,4]),as.numeric(power_5A_escape[2,4]))))

#we had the power to detect an effect size of 28.8%

#What sample size would we have needed with our variation to see an effect size of 10%


#use this function
effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

escape_5A_power_table <- sample_calc(power_5A_escape, effect_percent,pooled_var_5A_escape)

ggplot(subset(escape_5A_power_table, effect > 10), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Effect size (%)")+
  ylab("Sample size needed")+
  theme_bw()


#to calculate pooled variance
power_5A_CT3 <- subset(growth, event_short == "5A"|event_short =="CT3") %>% summarize(
  n = n(),
  var = var(V801)
)

#to calculate n
power_5A_CT3_2 <- subset(growth, event_short == "5A"|event_short =="CT3") %>% group_by(event_short) %>% summarize(
  n = n(),
  var = var(V801)
)

sqrt((21831723*2/33))*(qt(0.975,32))

#The difference would have had to have been 2343.035 cm3 in order to detect the difference between these events as significant
#we detected a difference of 2103.73 cm3
###############

write.csv(growth,file = "2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")

###table to summarize mean volume per event
Volume_summary <- data_frame()
event_volume <- as.data.frame(aggregate(V801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(V801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(V801 ~ event_short, data = growth, FUN = length)

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

write.csv(event_volume_summary, file = "LC_2023/2023_growth_inventory_analysis/event_volume_table.csv")
event_volume_summary$construct = c("LC-102","LC-102","LC-102","Escape","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","LC-102","Escape","Control")

myColors <- c('gray',
              'cyan4',
              'red4')


names(myColors) <- levels(event_volume_summary$construct)
colorScale <- scale_color_manual(name = "construct",values = myColors)


volume_plot <- ggplot(event_volume_summary, aes(x=reorder(Event,mean_Volume),y=mean_Volume, color = construct))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_volume,ymax=upper_volume))+
  xlab("Event")+
  ylab("Volume index (cm^3)")+
  theme_bw()

volume_plot_2 <- volume_plot + colorScale

ggsave(filename = "LC_2023/2023_growth_inventory_analysis/volume_plot_2023.png", plot = volume_plot_2, dpi = 300)

###Diameter summary
diameter_summary <- data_frame()
event_diameter <- as.data.frame(aggregate(D801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(D801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(D801 ~ event_short, data = growth, FUN = length)

event_diameter_summary <- inner_join(event_diameter, event_var, by = "event_short")
event_diameter_summary <- inner_join(event_diameter_summary, event_n, by = "event_short")

colnames(event_diameter_summary) <- c("Event", "mean_diameter","var","n")

event_diameter_summary$se <- sqrt(event_diameter_summary$var/event_diameter_summary$n)
event_diameter_summary$tst <- qt(0.975, event_diameter_summary$n)
event_diameter_summary$sd <- sqrt(event_diameter_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_diameter_summary$upper_diameter <- event_diameter_summary$mean_diameter + (event_diameter_summary$tst*event_diameter_summary$se)
event_diameter_summary$lower_diameter <- event_diameter_summary$mean_diameter - (event_diameter_summary$tst*event_diameter_summary$se)

write.csv(event_diameter_summary, file = "LC_2023/2023_growth_inventory_analysis/event_diameter_table.csv")

ggplot(event_diameter_summary, aes(x=reorder(Event,mean_diameter),y=mean_diameter))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_diameter,ymax=upper_diameter))+
  xlab("Event")+
  ylab("diameter (mm)")


##Height summary#####
height_summary <- data_frame()
event_height <- as.data.frame(aggregate(H801 ~ event_short, data = growth, FUN = mean))
event_var <- aggregate(H801 ~ event_short, data = growth, FUN = var)
event_n <- aggregate(H801 ~ event_short, data = growth, FUN = length)

event_height_summary <- inner_join(event_height, event_var, by = "event_short")
event_height_summary <- inner_join(event_height_summary, event_n, by = "event_short")

colnames(event_height_summary) <- c("Event", "mean_height","var","n")

event_height_summary$se <- sqrt(event_height_summary$var/event_height_summary$n)
event_height_summary$tst <- qt(0.975, event_height_summary$n)
event_height_summary$sd <- sqrt(event_height_summary$var)
qt(0.975, 21)
qt(0.975, 29)
#can estimate quartile for all to be 2
event_height_summary$upper_height <- event_height_summary$mean_height + (event_height_summary$tst*event_height_summary$se)
event_height_summary$lower_height <- event_height_summary$mean_height - (event_height_summary$tst*event_height_summary$se)

write.csv(event_height_summary, file = "LC_2023/2023_growth_inventory_analysis/event_height_table.csv")

ggplot(event_height_summary, aes(x=reorder(Event,mean_height),y=mean_height))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_height,ymax=upper_height))+
  xlab("Event")+
  ylab("height (mm)")





###Volume w/o transformation

# Fit model to test a difference between transgenic and control (construct2)
# Check block effect
# Volume
mod000 <- lm(V801 ~ V49 + construct2 + block, data = growth)
summary(mod000)
#Block is not significant but has influence

# Make a block as a random effect
v_mod3.1 <- lme(V801 ~ V49 + construct2, random = ~1|block, data = growth)
summary(v_mod3.1)

# Residual and qq plots
plot(fitted(v_mod3.1), residuals(v_mod3.1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod3.1)); qqline(residuals(v_mod3.1))

plot(v_mod3.1, which=1)



v_mod4.1 <- lme(V801 ~ V49 + event2, random = ~1|block, data=growth)
summary(v_mod4.1)

plot(fitted(v_mod4.1), residuals(v_mod4.1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(v_mod4.1)); qqline(residuals(v_mod4.1))
#definitely don't look as good as height and diameter.
#Q-Q plot shows there is a prominent right tail.

emmeans(v_mod4.1, specs = pairwise ~event2)

predict_set_v_III <- subset(growth, select = c("ID","event","event_short","construct","construct2","block","V49","V801"))
predict_set_v_III$predicted_V801 <- predict(v_mod4.1, new_data = predict_set_v_III, interval = "confidence")

#check
ggplot(predict_set_v_III,aes(x= V801, y = predicted_V801))+
  geom_point()
mean(predict_set_v_III$V801)

#average modeled volume
modelled_CT3_mean_v_II <- mean(predict_set_v_III$predicted_V801)
modelled_pooled_escape_mean_v_II <- mean(predict_set_v_III$predicted_V801) + 1942923
modelled_5A_mean_v_II <- mean(predict_set_v_III$predicted_V801) + 3973925
modelled_5C_mean_v_II <- mean(predict_set_v_III$predicted_V801) + 2073592
modelled_5_mean_v_II <- mean(predict_set_v_III$predicted_V801) + 1202274
modelled_2H_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 244672
modelled_13_15E_mean_v_II <- mean(predict_set_v_III$predicted_V801) + 714406
modelled_13_15B_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 1079120
modelled_4A_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 2547238 
modelled_4B_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 1550164
modelled_1_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 996201
modelled_1C_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 1129050
modelled_7_mean_v_II <- mean(predict_set_v_II$predicted_V801) + 1567743

event_volume_effect_II <- as.data.frame(c(modelled_CT3_mean_v_II, modelled_pooled_escape_mean_v_II, modelled_5A_mean_v_II, modelled_5C_mean_v_II, modelled_5_mean_v_II, modelled_2H_mean_v_II, modelled_13_15E_mean_v_II, modelled_13_15B_mean_v_II, modelled_4A_mean_v_II, modelled_4B_mean_v_II, modelled_1_mean_v_II, modelled_1C_mean_v_II, modelled_7_mean_v_II))

row.names(event_volume_effect_II) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B","4A","4B","1","1C","7")
colnames(event_volume_effect_II) <- c("volume")


###############

