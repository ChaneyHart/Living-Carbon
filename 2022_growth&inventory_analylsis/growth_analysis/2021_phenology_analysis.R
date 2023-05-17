#Budflush phenology analysis from Spring, 2022

#open packages (may need to install them first)
library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)

#set working directory to whereever raw data is located. This will be different on your computer.
setwd("~/Google Drive/G_Living Carbon")
#import data 
LB_budflush_all <- read.csv("LC_2022_Jan-May/Budflush/LB_budflush_all3.csv")

#clean up
phen <- na.omit(LB_budflush_all)
#phen <- subset(LB_budflush_all[,1:34])
phen$X7 <- as.numeric(phen$X7)
phen$X10 <- as.numeric(phen$X10)
phen$X14 <- as.numeric(phen$X14)
phen$X16 <- as.numeric(phen$X16)
phen$X21 <- as.numeric(phen$X21)

#######
#convert datatable to a longer format.
#columns for each week will collapsed and their sums will be added into a new column called "Day"

#phen_long <- pivot_longer(phen, cols = c(day0,day7,day10,day14,day16,day21), names_to = "Day", values_to = "Score")
phen_long <- pivot_longer(phen, cols = c(X0,X7,X10,X14,X16,X21), names_to = "Day", values_to = "Score")

#Top.long <- pivot_longer(phen, cols = c(Top_0,Top_7,Top_10,Top_14,Top_16,Top_21), names_to = "Day", values_to = "Top_score")
#Middle.long <- pivot_longer(phen, cols = c(Middle_0,Middle_7,Middle._10,Middle_14,Middle_16,Middle_21), names_to = "Day", values_to = "Middle_score")
#Bottom.long <- pivot_longer(phen, cols = c(Bottom_0,Bottom_7,Bottom_10,Bottom_14,Bottom_16,Bottom_21), names_to = "Day", values_to = "Bottom_score")
#Avg.long <- pivot_longer(phen, cols = c(Avg_0,Avg_7,Avg_10,Avg_14,Avg_16,Avg_21), names_to = "Day", values_to = "Avg_score")



#function to group data by individual plant
by_plant <- group_by(sen_long, ID)
day2<- substr(by_plant$Day,2,3)
by_plant$day2 <- as.numeric(day2)


#by_plant$day2 <- as.numeric(by_plants$Day)
#by_plant$day2 <- 
#by_plant <- group_by(Top.long, ID)
#by_plant2 <- group_by(Middle.long, ID)
#by_plant3 <- group_by(Bottom.long, ID)
#by_plant4 <- group_by(Avg.long, ID)

#nested1 <- by_plant %>%
# nest(data = ID)
#plot out data
plot(Score~day2, data= by_plant)

#fit linear model (regression line) for each plant's bud scores.
#allows to estimate at which day buds reached a certain development point for each tree

models <- do(by_plant, 
             bind_rows(coef(lm(Score~day2, data = .))))


#fit lm for each plant
#models <- do(by_plant,
#  bind_rows(coef(lm(Top_score~Day, data = .))))
#models2 <- do(by_plant2,
#  bind_rows(coef(lm(bottom~day, data = .))))
#models3 <- do(by_plant3,
# bind_rows(coef(lm(bottom~day, data = .))))
#models4 <- do(by_plant4,
#             bind_rows(coef(lm(bottom~day, data = .))))

#phen.top <- subset(LB_budflush_all, select = c(ID,construct,event,Top_0,Top_7,Top_10,Top_14,Top_16,Top_21))
#phen.middle <- subset(LB_budflush_all, select = c(ID,construct,event,Middle_0,Middle_7,Middle._10,Middle_14,Middle_16,Middle_21))
#phen.bottom <- subset(LB_budflush_all, select = c(ID,construct,event,Bottom_0,Bottom_7,Bottom_10,Bottom_14,Bottom_16,Bottom_21))
#phen.avg <- subset(LB_budflush_all, select = c(ID,construct,event,Avg_0,Avg_7,Avg_10,Avg_14,Avg_16,Avg_21))

#quickplot(ID, Top_14, data = phen.top, color = construct)
#quickplot(ID, Middle_14, data = phen.middle, color = construct)
#quickplot(ID, Bottom_14, data = phen.bottom, color = construct)
#quickplot(ID, Top_7, data = phen.top, color = construct)
#quickplot(ID, Top_10, data = phen.top, color = construct)
#quickplot(ID, Top_16, data = phen.top, color = construct)

# set rating threshold
#This is the score associated with the development phase of interest
rating <- 2.5

#this is y = mx+b rearranged to solve for x (days) when y = rating
models$days_est <- (rating - models$`(Intercept)`)/models$day2 
#plot(Score~day2, data=by_plant$ID`$"LCOR-001"`)
#models2$days_est <- (rating - models2$`(Intercept)`)/models2$day2 

#convert to factor
models$id <- as.factor(models$ID)

#sort data by days_est
models <- models[order(models$days_est),]
#remove negative values
models <- models[(-1:-6),]
#models2 <- models2[order(models2$days_est),]

#create index column
models$index <- as.numeric(row.names(models))
#models2$index <- as.numeric(row.names(models2))

write.table(models, "~/phenology_models_average.txt", sep="\t") 
#write.table(models2, "C:/Users/amand/Desktop/senbottom.txt", sep="\t") 


#plot out data
ggplot(models, aes(x=index, y=days_est, color=id)) +
  geom_point()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(models2, aes(x=index, y=days_est, color=id)) +
  geom_point()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#######some extra stuff######

#create new df with model statistics for each plant
model_stats <- do(by_plant,
                  glance(lm(Score ~ day2, data = .)))

#make new df using information for models to augment each observation (predictions, residuals, standard errors)
augmented <- do(by_plant,
                augment(
                  lm(Score~day2, data=.)
                ))

#plot out residuals versus day
ggplot(augmented, aes(day2, .resid)) +
  geom_line(aes(group = ID), alpha = 1 / 3) + 
  geom_smooth(se = FALSE) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
