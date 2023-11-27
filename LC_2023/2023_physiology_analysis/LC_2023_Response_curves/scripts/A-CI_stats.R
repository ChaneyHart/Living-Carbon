#stats A-Ci parameters
# Additional libraries
library(emmeans)
library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)
#import cleaned ACI data

ACI_summary_tree <- read.csv("LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_summary_tree_cleaned.csv")
ACI_summary_tree$PhiPS2 <- as.numeric(ACI_summary_tree$PhiPS2)
ACI_summary_tree$ETR <- as.numeric(ACI_summary_tree$ETR)

#import photosynthetic parameters
ACI_parameters <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/aci_parameters_list.csv")

#stats

unique(aci_parameters_list$event)

# Create tiers factor
ACI_parameters <- ACI_parameters %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~ "Control",
  event_short == "CT3" ~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))

##########################
#####Vcmax
# 

ACI_parameters$sample_date <- as.factor(ACI_parameters$sample_date)
# summary stats
mean_vcmax <- mean(ACI_parameters$Vcmax)
std_Vcmax <- sd(ACI_parameters$Vcmax)
CV_vcmax <- (std_Vcmax/mean_vcmax)*100

#fit model
mod0 <- lm(Vcmax ~ tier + sample_date + block, data = ACI_parameters)
#fit aov
aov0 <- aov(Vcmax ~ tier + sample_date + block, data = ACI_parameters)
summary(mod0)
anova(aov0)

# Make a sample_date as a random effect
mod1 <- lme(Vcmax ~ tier, random = ~1|sample_date, data = ACI_parameters)
summary(mod1)

# Residual and qq plots
plot(fitted(mod1), residuals(mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod1)); qqline(residuals(mod1))
#one strong outlier: residual = -52.9
fitted(mod1)
residuals(mod1)

#pairwise comparison


emmeans_Vcmax <- emmeans(mod1, specs = pairwise ~tier)
emmeans_Vcmax

# graph modelled tier means
tier_Vcmax_effect <- as.data.frame(emmeans_Vcmax$emmeans)
tier_Vcmax_effect$modelled_Vcmax <- (tier_Vcmax_effect$emmean)

my_colors3 <- c("gray0","chartreuse4","indianred3","gray")
modelled_vcmax_plot <- ggplot(tier_Vcmax_effect, aes(x=tier,y=emmean,color=tier))+
  geom_point(size =2)+
  ylab("Modelled Vcmax (µmol/m2/sec)")+
  geom_errorbar(aes(ymin=emmean-2*SE,ymax=emmean+2*SE),linewidth=1)+
  scale_color_manual(values = my_colors3)+
  theme_bw()

modelled_vcmax_plot

ggsave(filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/modelled_Vcmax.png",plot= modelled_vcmax_plot, dpi=300)

####what effect size would have been needed to observe sig results when comparing top events to escapes


power_top_control_Vcmax <- subset(ACI_parameters, tier == "top"|tier =="Control") %>% group_by(tier) %>% summarize(
  n = n(),
  var = var(Vcmax),
  mean = mean(Vcmax),
  sd = sd(Vcmax),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

top_control_Vcmax_diff_mean <- (power_top_control_Vcmax[2,4] - power_top_control_Vcmax[1,4])
top_control_Vcmax_diff_lower <- (power_top_control_Vcmax[2,8] - power_top_control_Vcmax[1,8])
top_control_Vcmax_diff_upper <- (power_top_control_Vcmax[2,9] - power_top_control_Vcmax[1,9])
pooled_var_top_control_Vcmax <- ((power_top_control_Vcmax[1,2]*power_top_control_Vcmax[1,3]) + (power_top_control_Vcmax[2,2]*power_top_control_Vcmax[2,3]))/(power_top_control_Vcmax[1,2]+power_top_control_Vcmax[2,2]-2)

qt(0.975,8)
sig_effect_size_top_control_Vcmax <- sqrt(((pooled_var_top_control_Vcmax[1,1])*2)/8)*(qt(0.975,8))/(mean(c(as.numeric(power_top_control_Vcmax[1,4]),as.numeric(power_top_control_Vcmax[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
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


effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

top_control_Vcmax_power_table <- sample_calc(power_top_control_Vcmax, effect_percent,pooled_var_top_control_Vcmax)

Vcmax_power_graph <- ggplot(subset(top_control_Vcmax_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=8), color="blue",linetype = 2)+
  geom_line(aes(x=21.2), color = "red", linetype = 2)

Vcmax_power_graph
ggsave(Vcmax_power_graph, filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/Vcmax_power_graph.png",dpi = 300)


#####
#Jmax

mean_Jmax <- mean(ACI_parameters$Jmax)
std_Jmax <- sd(ACI_parameters$Jmax)
CV_Jmax <- (std_Jmax/mean_Jmax)*100

#fit model
mod00 <- lm(Jmax ~ tier + sample_date + block, data = ACI_parameters)
#fit aov
aov00 <- aov(Jmax ~ tier + sample_date + block, data = ACI_parameters)
summary(mod00)
anova(aov00)

# Make a sample_date as a random effect
mod2 <- lme(Jmax ~ tier, random = ~1|sample_date, data = ACI_parameters)
summary(mod2)

# Residual and qq plots
plot(fitted(mod2), residuals(mod2), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod2)); qqline(residuals(mod2))
#one strong outlier: residual = -52.9
fitted(mod2)
residuals(mod2)

#pairwise comparison


emmeans_Jmax <- emmeans(mod2, specs = pairwise ~tier)
emmeans_Jmax

# graph modelled tier means
tier_Jmax_effect <- as.data.frame(emmeans_Jmax$emmeans)
tier_Jmax_effect$modelled_Jmax <- (tier_Jmax_effect$emmean)

my_colors3 <- c("gray0","chartreuse4","indianred3","gray")
modelled_Jmax_plot <- ggplot(tier_Jmax_effect, aes(x=tier,y=emmean,color=tier))+
  geom_point(size =2)+
  ylab("Modelled Jmax (µmol/m2/sec)")+
  geom_errorbar(aes(ymin=emmean-2*SE,ymax=emmean+2*SE),linewidth=1)+
  scale_color_manual(values = my_colors3)+
  theme_bw()

modelled_Jmax_plot

ggsave(filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/modelled_Jmax.png",plot= modelled_Jmax_plot, dpi=300)

####what effect size would have been needed to observe sig results when comparing top events to escapes


power_top_control_Jmax <- subset(ACI_parameters, tier == "top"|tier =="Control") %>% group_by(tier) %>% summarize(
  n = n(),
  var = var(Jmax),
  mean = mean(Jmax),
  sd = sd(Jmax),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

top_control_Jmax_diff_mean <- (power_top_control_Jmax[2,4] - power_top_control_Jmax[1,4])
top_control_Jmax_diff_lower <- (power_top_control_Jmax[2,8] - power_top_control_Jmax[1,8])
top_control_Jmax_diff_upper <- (power_top_control_Jmax[2,9] - power_top_control_Jmax[1,9])
pooled_var_top_control_Jmax <- ((power_top_control_Jmax[1,2]*power_top_control_Jmax[1,3]) + (power_top_control_Jmax[2,2]*power_top_control_Jmax[2,3]))/(power_top_control_Jmax[1,2]+power_top_control_Jmax[2,2]-2)

qt(0.975,8)
sig_effect_size_top_control_Jmax <- sqrt(((pooled_var_top_control_Jmax[1,1])*2)/8)*(qt(0.975,8))/(mean(c(as.numeric(power_top_control_Jmax[1,4]),as.numeric(power_top_control_Jmax[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
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


effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

top_control_Jmax_power_table <- sample_calc(power_top_control_Jmax, effect_percent,pooled_var_top_control_Jmax)

Jmax_power_graph <- ggplot(subset(top_control_Jmax_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=8), color="blue",linetype = 2)+
  geom_line(aes(x=13), color = "red", linetype = 2)

Jmax_power_graph
ggsave(Jmax_power_graph, filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/Jmax_power_graph.png",dpi = 300)

######
#Amax

mean_Amax <- mean(ACI_parameters$Amax)
std_Amax <- sd(ACI_parameters$Amax)
CV_Amax <- (std_Amax/mean_Amax)*100

#fit model
mod0000 <- lm(Amax ~ tier + sample_date + block, data = ACI_parameters)
#fit aov
aov0000 <- aov(Amax ~ tier + sample_date + block, data = ACI_parameters)
summary(mod0000)
anova(aov0000)

# Make a sample_date as a random effect
mod3 <- lme(Amax ~ tier, random = ~1|sample_date, data = ACI_parameters)
summary(mod3)

# Residual and qq plots
plot(fitted(mod3), residuals(mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod3)); qqline(residuals(mod3))


#pairwise comparison


emmeans_Amax <- emmeans(mod3, specs = pairwise ~tier)
emmeans_Amax

# graph modelled tier means
tier_Amax_effect <- as.data.frame(emmeans_Amax$emmeans)
tier_Amax_effect$modelled_Amax <- (tier_Amax_effect$emmean)

my_colors3 <- c("gray0","chartreuse4","indianred3","gray")
modelled_Amax_plot <- ggplot(tier_Amax_effect, aes(x=tier,y=emmean,color=tier))+
  geom_point(size =2)+
  ylab("Modelled Amax (µmol/m2/sec)")+
  geom_errorbar(aes(ymin=emmean-2*SE,ymax=emmean+2*SE),linewidth=1)+
  scale_color_manual(values = my_colors3)+
  theme_bw()

modelled_Amax_plot

ggsave(filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/modelled_Amax.png",plot= modelled_Amax_plot, dpi=300)

####what effect size would have been needed to observe sig results when comparing top events to escapes


power_top_control_Amax <- subset(ACI_parameters, tier == "top"|tier =="Control") %>% group_by(tier) %>% summarize(
  n = n(),
  var = var(Amax),
  mean = mean(Amax),
  sd = sd(Amax),
  se = sd/sqrt(n),
  tst = qt(0.975, n-1),
  lower = mean - tst*se,
  upper = mean + tst*se
)

top_control_Amax_diff_mean <- (power_top_control_Amax[2,4] - power_top_control_Amax[1,4])
top_control_Amax_diff_lower <- (power_top_control_Amax[2,8] - power_top_control_Amax[1,8])
top_control_Amax_diff_upper <- (power_top_control_Amax[2,9] - power_top_control_Amax[1,9])
pooled_var_top_control_Amax <- ((power_top_control_Amax[1,2]*power_top_control_Amax[1,3]) + (power_top_control_Amax[2,2]*power_top_control_Amax[2,3]))/(power_top_control_Amax[1,2]+power_top_control_Amax[2,2]-2)

qt(0.975,8)
sig_effect_size_top_control_Amax <- sqrt(((pooled_var_top_control_Amax[1,1])*2)/8)*(qt(0.975,8))/(mean(c(as.numeric(power_top_control_Amax[1,4]),as.numeric(power_top_control_Amax[2,4]))))

#What sample size would we have needed with our variation to see an effect size of 10%

#use this function
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


effect_percent <- as.data.frame(c(5,10,15,20,25,30,35,40,45,50))

top_control_Amax_power_table <- sample_calc(power_top_control_Amax, effect_percent,pooled_var_top_control_Amax)

Amax_power_graph <- ggplot(subset(top_control_Amax_power_table, effect > 5), aes(x=effect, y = n))+
  geom_point()+
  geom_line()+
  xlab("Real population effect size (µ2-µ1)/µ1 (%)")+
  ylab("Sample size needed")+
  theme_bw()+
  geom_line(aes(y=8), color="blue",linetype = 2)+
  geom_line(aes(x=12), color = "red", linetype = 2)

Amax_power_graph
ggsave(Amax_power_graph, filename = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/Amax_power_graph.png",dpi = 300)


