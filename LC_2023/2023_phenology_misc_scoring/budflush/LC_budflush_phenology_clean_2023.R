#Budflush phenology analysis

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)

#setwd("C:/Users/Veron/OneDrive/Documents/URSA")
#import
LB_budflush_all <- read.csv("LC_2023/2023_phenology_misc_scoring/2023_budflush_all.csv")
LB_budflush_all <- na.omit(LB_budflush_all)
LB_budflush_all <- mutate(LB_budflush_all, Avg_0 = rowMeans(select(LB_budflush_all, c(Top_0, Bottom_0)), na.rm = TRUE))
LB_budflush_all <- mutate(LB_budflush_all, Avg_2 = rowMeans(select(LB_budflush_all, c(Top_2, Bottom_2)), na.rm = TRUE))
LB_budflush_all <- mutate(LB_budflush_all, Avg_6 = rowMeans(select(LB_budflush_all, c(Top_6, Bottom_6)), na.rm = TRUE))

LB_budflush_all
#############################

#subset into top, bottom, average
phen.top <- subset(LB_budflush_all, select = c(ID,Top_0,Top_2,Top_6))
phen.bottom <- subset(LB_budflush_all, select = c(ID,Bottom_0,Bottom_2,Bottom_6))
phen.avg <- subset(LB_budflush_all, select = c(ID,Avg_0,Avg_2,Avg_6))
#clean up
phen.top <- na.omit(phen.top)
phen.bottom <- na.omit(phen.bottom)
phen.avg <- na.omit(phen.avg)
#check data types
str(phen.top)
str(phen.bottom)
str(phen.avg)
#change to numeric
phen.avg$Avg_0 <- as.numeric(phen.avg$Avg_0)
phen.avg$Avg_2 <- as.numeric(phen.avg$Avg_2)
phen.avg$Avg_6 <- as.numeric(phen.avg$Avg_6)

phen.avg <- na.omit(phen.avg)
#######

#convert to long format
phen.top_long <- pivot_longer(phen.top, cols = c(Top_0,Top_2,Top_6), names_to = "Day", values_to = "Score")
#phen.top_long <- phen.top_long %>% select(ID,event,construct,Day,Score) %>% mutate(day2 = as.numeric(substr(phen.top_long$Day,5,7)))
phen.top_long <- phen.top_long %>% mutate(day2 = as.numeric(substr(phen.top_long$Day,5,7)))

phen.bottom_long <- pivot_longer(phen.bottom, cols = c(Bottom_0,Bottom_2,Bottom_6), names_to = "Day", values_to = "Score")
phen.bottom_long <- phen.bottom_long  %>% mutate(day2 = as.numeric(substr(phen.bottom_long$Day,8,10)))

phen.avg_long <- pivot_longer(phen.avg, cols = c(Avg_0,Avg_2,Avg_6), names_to = "Day", values_to = "Score")
phen.avg_long <- phen.avg_long %>%  mutate(day2 = as.numeric(substr(phen.avg_long$Day,5,7)))


#function to group by plant
by_plant <- group_by(phen.top_long, ID)

by_plant3 <- group_by(phen.bottom_long, ID)
by_plant4 <- group_by(phen.avg_long, ID)

#fit lm for each plant
models_top <- do(by_plant, 
             bind_rows(coef(lm(Score~day2, data = .))))


models_bottom <- do(by_plant3,
  bind_rows(coef(lm(Score~day2, data = .))))
models_avg <- do(by_plant4,
  bind_rows(coef(lm(Score~day2, data = .))))

# set rating threshold
rating <- 2.5

#this is y = mx+b rearranged to solve for x (days) when y = rating
models_top$days_est <- (rating - models_top$`(Intercept)`)/models_top$day2 
 
models_bottom$days_est <- (rating - models_bottom$`(Intercept)`)/models_bottom$day2 
models_avg$days_est <- (rating - models_avg$`(Intercept)`)/models_avg$day2 

#convert to factor
models_top$ID <- as.factor(models_top$ID)

models_bottom$ID <- as.factor(models_bottom$ID)
models_avg$ID <- as.factor(models_avg$ID)

#sort data by days_est
models_top <- models_top[order(models_top$days_est),]

models_bottom <- models_bottom[order(models_bottom$days_est),]
models_avg <- models_avg[order(models_avg$days_est),]

#set negative values to 0

models_top$days_est[models_top$days_est<0] <- 0
models_bottom$days_est[models_bottom$days_est<0] <- 0
models_avg$days_est[models_avg$days_est<0] <- 0


#create index column
models_top$index <- as.numeric(row.names(models_top))
models_bottom$index <- as.numeric(row.names(models_bottom))
models_avg$index <- as.numeric(row.names(models_avg))

write.table(models_top, "~/phenology_models_top_2023.txt", sep="\t") 
write.table(models_bottom, "~/phenology_models_bottom_2023.txt", sep="\t") 
write.table(models_avg, "~/phenology_models_average_2023.txt", sep="\t") 


#plot out data
ggplot(models_top, aes(x=index, y=days_est, color=ID)) +
  geom_point()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(models_bottom, aes(x=index, y=days_est, color=ID)) +
  geom_point()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

phen_dist <- ggplot(models_avg, aes(x=index, y=days_est, color=ID)) +
  geom_point()+
  xlab("tree")+
  ylab("days till score = 2.5")+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

phen_dist
ggsave("LC_2023/phen_dist_2023.png", plot = phen_dist, dpi = 300)

########## statistical comparison between events###########

#bring in growth data w/ info about each tree

growth <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")



growth_phen_compare_inner <- inner_join(growth, models_avg, by = "ID")




###phenology analysis#####
phenology_comp <- ggplot(growth_phen_compare_inner, aes(x = reorder(event_short,days_est), y = days_est, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("avg days till budflush = 2.5")


phenology_comp
ggsave("bud.phenology_2023.png", plot =phenology_comp, width = 10, height = 5, dpi =300, units = "in")


anova_phenology <- lm(days_est ~ event, data = growth_phen_compare_inner)
aov_phenology <- aov(days_est ~ event, data = growth_phen_compare_inner)
#plot residuals
plot(aov_phenology,which = 1)
# datapoint 394 (LCOR-555) looks like an outlier.
# if the standardized residuals of a datapoint are greater than 2.0 or less than -2.0, that point is recommended to be removed.
plot(aov_phenology,which = 3)

#standardized residual for point 394 is greater than 2, remove from dataset and redefine outliers

#model without point 394
growth_phen_compare_inner_II <- growth_phen_compare_inner[-394,]
anova_phenology_II <- lm(days_est ~ event, data = growth_phen_compare_inner_II)
aov_phenology_II <- aov(days_est ~ event, data = growth_phen_compare_inner_II)
plot(aov_phenology_II,which = 1)
plot(aov_phenology_II,which = 3)

plot(aov_phenology_II,which = 2)
plot(aov_phenology_II,which = 4)


##Stats comparison between construct and events using the model just formed

# Additional libraries to load (may need to install packages)
library(emmeans)
library(nlme)
library(lme4)

# Combine 16-20 and 8-9D (This will be useful later)
unique(growth_phen_compare_inner_II$event)


growth_phen_compare_inner_II$event2 <- growth_phen_compare_inner_II$event
growth_phen_compare_inner_II$event2[growth_phen_compare_inner_II$event == "LC-102 16-20" |
                growth_phen_compare_inner_II$event == "LC-102 8-9D"] <- "escape"

######Construct model(s)#########

# Fit model to test a difference between transgenic and control (construct2)
# Check block effect
# Height
phen_mod0 <- lm(days_est ~ construct2 + block, data = growth_phen_compare_inner_II)
summary(phen_mod0)
#block is a significant predictor of phenology but not of primary interest. Event is of interest. Therefore, the block effect will be treated as a random effect.

# Make a block as a random effect
phen_mod1 <- lme(days_est ~ construct2, random = ~1|block, data = growth_phen_compare_inner_II)
summary(phen_mod1)

# Residual and qq plots to check assumptions of model
# 1 - residuals test that different groups have equal variance
plot(fitted(phen_mod1), residuals(phen_mod1), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
#Variance of third grouping is more spread, but assumption is not egregiously violated
qqnorm(residuals(phen_mod1)); qqline(residuals(phen_mod1))
#This does not look greaat, let's look at something different
##

#ok, in this model we are looking at differences between transgenic, control and escape
phen_mod3 <- lme(days_est ~ construct, random = ~1|block, data = growth_phen_compare_inner_II)
summary(phen_mod3)

# Residual and qq plots
plot(fitted(phen_mod3), residuals(phen_mod3), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(phen_mod3)); qqline(residuals(phen_mod3))
#still don't love how this looks but let's look at results of this model.

#This table documents estimates of differences between contstruct groups along with their standard error (and associated p-values)
emmeans(phen_mod3, specs = pairwise ~construct)
#I interpret this as: Difference between control and escape is minimal and not statistically significant (not surprising)
#Difference between LC-102 (transgenic) and both control and escape is largish (about a day) and is statistical significant.
#While there was perhaps some non-normality in the data, this result was pretty strong (not a borderline call) - so I tend to believe it.


######### Effect of event (escapes pooled) ###########
#again, using block as random effect
phen_mod4 <- lme(days_est ~ event2, random = ~1|block, data = growth_phen_compare_inner_II)
summary(phen_mod4)

# Residual and qq plots to check assumptions again
plot(fitted(phen_mod4), residuals(phen_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)

#This looks pretty good, most groups have equal variance (spread on either side of 0).
qqnorm(residuals(phen_mod4)); qqline(residuals(phen_mod4))
#This looks much better than the case for the last model. The data points fall closer to the line.

#Again statistical comparison table is generated. This time the difference between events is estimated from the model as well as standard error and associated p-value.
emmeans(phen_mod4, specs = pairwise ~event2)
#There are many comparisons that are significant! Interesting!
#These can be dug into more.


##overal conclusions###
#There is evidence of a difference in the date of budflush between various transgenic events and various control events.
#In addition, there is evidence that as a whole there is a difference between all transgenic events and all control events and between all transgenic events and all escape events.



#######Effect of phenology on heat damage in May heat wave#######

May_heat_scores <- read.csv("LC_2023/May_heat_scoring.csv", skip = 13)

Heat_score <- inner_join(growth_phen_compare_inner, May_heat_scores, by = "ID")


ggplot(Heat_score, aes(x=days_est, y = Score, color = event_short))+
  geom_point()

#I thought this may be interesting - don't see much of a trend.
