#Biomass monte carlo

library(ggplot2)
library(dplyr)
library(readxl)
#install.packages("fabricatr")
library(fabricatr)
library(purrr)
library(tidyr)
library(nlme)


growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
#growth_dat$Vol_index <- 0.333*((3.14*((growth_dat$DBH/2)^2))*growth_dat$Height)

growth_dat <- growth_dat %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~"control", 
  event_short == "CT3" ~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

growth_dat <- growth_dat %>% mutate(section = case_when(
  column <= 8 & row <= 8 ~ 1,
  column <= 8 & row > 8 ~ 2,
  column > 8 & column <= 16 & row <= 8 ~ 3,
  column > 8 & column <= 16 & row > 8 ~ 4,
  column > 16 & column <= 24 & row <= 8 ~5,
  column > 16 & column <= 24 & row > 8 ~ 6,
  column > 24 & column <= 28 & row <= 8 ~7,
  column > 24 & column <= 28 & row > 8 ~8,
  column > 44 & column <= 50 & row <=8 ~9,
  column > 44 & column <= 50 & row > 8 ~10,
  column > 50 & column < 56 & row <= 8 ~ 11,
  column > 50 & column < 56 & row > 8 ~ 12
))


#Split each event into terciles

biomass_set <- growth_dat %>% group_by(event_short) %>% mutate(
  strata = ntile(V801,3))

#Select only events of interest
biomass_subset <- subset(biomass_set, event_short == "13-15E" | event_short == "5A" | event_short == "4A" | event_short == "16-20" | event_short == "8-9D" | event_short == "CT3")

# add info about the total number of trees in each event (N)
biomass_subset <- biomass_subset %>% mutate(N= case_when(
  event_short == "13-15E"~"32",
  event_short == "16-20"~"33",
  event_short == "5A"~"34",
  event_short == "4A"~"19",
  event_short == "2H"~ "29",
  event_short == "5C"~ "33",
  event_short == "8-9D"~"29",
  event_short == "CT3"~"34"
))



biomass_subset$N <- as.numeric(biomass_subset$N)
#for each event, we want 50% of the tree samples
biomass_subset$n_final <- biomass_subset$N*0.5

#testing function

biomass_5A <- subset(biomass_subset, event_short == "5A")
biomass_4A <- subset(biomass_subset, event_short == "4A")
biomass_13_15E <- subset(biomass_subset, event_short == "13-15E")
biomass_16_20 <- subset(biomass_subset, event_short == "16-20")
biomass_8_9D <- subset(biomass_subset, event_short == "8-9D")
biomass_CT3 <- subset(biomass_subset, event_short == "CT3")

#dat_5A <- subset(biomass_5A, select = c("ID","strata","V801","denom","N","n_final"))
#colnames <- c("ID","strata", "V801","denom","N","n_final")
#colnames <- match(colnames, names(dat_5A))
#ID <- (as.matrix(dat_5A[,colnames[1]]))
#Strata <- (as.matrix((dat_5A[,colnames[2]])))
#Volume <- (as.matrix(dat_5A[,colnames[3]]))
#denom <- (as.matrix(dat_5A[,colnames[4]]))
#N <- (as.matrix(dat_5A[1,colnames[5]]))
#n <- (as.matrix(dat_5A[1,colnames[6]]))

#SV <- mean(Volume)
#SD = sd(Volume)

#V <- matrix(nrow=1000, ncol=5)
#B <- double(5) 
#MSE <- double(5)

## SRS
#samp <- sample(Volume, n) 
#V[1,1] <- sum(samp/(n/N))/N
#B[1] <- B[1] + (V[1,1] - SV)
#MSE[1] <- MSE[1] + (V[1,1] - SV)^2

## STRS

#strt <- levels(factor(Strata))
#Nh <- gapply(data.frame(Strata), FUN=function(x) nrow(x), form=~Strata) #Length of each strata
#nh <- c(2,2,11) #sample size of each strata with neyman allocation

#V[1,2] <- 0
#sample(1:Nh[3],nh[3])
#Volume[Strata==strt[1]][c(7,2,4)]

#for(j in 1:length(Nh)){
 # samp <- sample(1:Nh[j], nh[j])
  #V[1,2] <- (V[1,2] + sum(Volume[Strata==strt[j]][samp])/(nh[j]/Nh[j]))/N
#}
#B[2] <- B[2] + (V[2,2] - SV)
#MSE[2] <- MSE[2] + (V[2,2] - SV)^2


##################################################################################

#with a dataset for a single event...
mcs_biomass <- function(data, colnames=c("ID","strata", "V801"), n,nh, reps=1000) {
  ## Match up column names
  colnames <- match(colnames, names(data))
  ID <-  as.matrix(data[,colnames[1]])	## ID 
  strata <- as.matrix(data[,colnames[2]])	## Strata
  Volume <- as.matrix(data[,colnames[3]])	## Volume
  
  
  ## Calculate population values
  N = length(Volume)
  SV <- mean(Volume)   ## mean biomass
  SD <- sd(Volume)
  
  ## Create initial values
  V <- matrix(nrow=reps, ncol=2)
  B <- double(2) 
  MSE <- double(2)
  
  ## Run Monte Carlo simulation
  for(i in 1:reps)
  { 
    ## SRS
    samp <- sample(Volume, n) 
    V[i,1] <- (sum(samp/(n/N)))/N
    B[1] <- B[1] + (V[i,1] - SV)
    MSE[1] <- MSE[1] + (V[i,1] - SV)^2
    
    ## STRS (neymans)
    strt <- levels(factor(strata))
    Nh <- gapply(data.frame(strata), FUN=function(x) nrow(x), form=~strata) #Population total per strata
    
    V[i,2] <- 0
    for(j in 1:length(Nh))
    {
      samp <- sample(1:Nh[j], nh[j])
      V[i,2] <- V[i,2] + (sum((Volume[strata==strt[j]][samp])/(nh[j]/Nh[j])))/N
    }
    B[2] <- B[2] + (V[i,2] - SV)
    MSE[2] <- MSE[2] + (V[i,2] - SV)^2
  }
  B <- B/reps
  MSE <- MSE/reps
  RMSE <- sqrt(MSE)
    
  #create return list
  lst <- list(n=n, reps=reps,Bias=B, MSE=MSE, RMSE=RMSE, Vhat=V, Vobs=SV)
  return(lst)
  
}

mc_5A <- mcs_biomass(biomass_5A,n =17,nh=c(3,3,11),reps = 1000)
mc_4A <- mcs_biomass(biomass_4A,n =17,nh=c(2,2,6),reps = 1000)
mc_16_20 <- mcs_biomass(biomass_16_20,n =17,nh=c(2,3,11),reps = 1000)
mc_8_9D <- mcs_biomass(biomass_8_9D,n =17,nh=c(4,2,8),reps = 1000)
mc_CT3 <- mcs_biomass(biomass_CT3,n =17,nh=c(4,4,9),reps = 1000)
mc_13_15E <- mcs_biomass(biomass_13_15E,n =17,nh=c(3,3,10),reps = 1000)


event <- c("5A","4A","16-20","8-9D","CT3","13-15E")
bias_SRS <- as.matrix(list(mc_5A[["Bias"]][1],mc_4A[["Bias"]][1],mc_16_20[["Bias"]][1],mc_8_9D[["Bias"]][1],mc_CT3[["Bias"]][1],mc_13_15E[["Bias"]][1]))
bias_STRS <- as.matrix(list(mc_5A[["Bias"]][2],mc_4A[["Bias"]][2],mc_16_20[["Bias"]][2],mc_8_9D[["Bias"]][2],mc_CT3[["Bias"]][2],mc_13_15E[["Bias"]][2]))
MSE_SRS <- as.matrix(list(mc_5A[["MSE"]][1],mc_4A[["MSE"]][1],mc_16_20[["MSE"]][1],mc_8_9D[["MSE"]][1],mc_CT3[["MSE"]][1],mc_13_15E[["MSE"]][1]))
MSE_STRS <- as.matrix(list(mc_5A[["MSE"]][2],mc_4A[["MSE"]][2],mc_16_20[["MSE"]][2],mc_8_9D[["MSE"]][2],mc_CT3[["MSE"]][2],mc_13_15E[["MSE"]][2]))
RMSE_SRS <- as.matrix(list(mc_5A[["RMSE"]][1],mc_4A[["RMSE"]][1],mc_16_20[["RMSE"]][1],mc_8_9D[["RMSE"]][1],mc_CT3[["RMSE"]][1],mc_13_15E[["RMSE"]][1]))
RMSE_STRS <- as.matrix(list(mc_5A[["RMSE"]][2],mc_4A[["RMSE"]][2],mc_16_20[["RMSE"]][2],mc_8_9D[["RMSE"]][2],mc_CT3[["RMSE"]][2],mc_13_15E[["RMSE"]][2]))


sample_comp <- as.data.frame(cbind(bias_SRS,bias_STRS,MSE_SRS,MSE_STRS,RMSE_SRS,RMSE_STRS))
sample_comp$event <- c("5A","4A","16-20","8-9D","CT3","13-15E")
colnames(sample_comp) <- c("bias_SRS","bias_STRS","MSE_SRS","MSE_STRS","RMSE_SRS","RMSE_STRS","event")

sample_comp_bias <- pivot_longer(sample_comp,cols = c("bias_SRS","bias_STRS"), names_to = "sampling_design" , values_to = c("bias"))
sample_comp_bias <- sample_comp_bias[,c(5:7)]
sample_comp_bias$bias <- as.matrix(unlist(sample_comp_bias$bias))
sample_comp_bias$sampling_design <- substr(sample_comp_bias$sampling_design,6,9)
str(sample_comp_bias)

sample_comp_MSE <- pivot_longer(sample_comp,cols = c("MSE_SRS","MSE_STRS"), names_to = "sampling_design" , values_to = c("MSE"))
sample_comp_MSE <- sample_comp_MSE[,c(5:7)]
sample_comp_MSE$MSE <- as.matrix(unlist(sample_comp_MSE$MSE))
sample_comp_MSE$sampling_design <- substr(sample_comp_MSE$sampling_design,5,9)

sample_comp_RMSE <- pivot_longer(sample_comp,cols = c("RMSE_SRS","RMSE_STRS"), names_to = "sampling_design" , values_to = c("RMSE"))
sample_comp_RMSE <- sample_comp_RMSE[,c(5:7)]
sample_comp_RMSE$RMSE <- as.matrix(unlist(sample_comp_RMSE$RMSE))
sample_comp_RMSE$sampling_design <- substr(sample_comp_RMSE$sampling_design,6,9)

?pivot_longer()
bias_plot <- ggplot(sample_comp_bias, aes(x=event,fill=sampling_design))+
  geom_col(aes(y=bias))+
  geom_col(aes(y=bias))
bias_plot


ggplot(sample_comp_bias, aes(event,abs(bias), fill = sampling_design))+
  geom_bar(stat = "identity",position = "dodge")+
  ylab("bias (absolute value)")

ggplot(sample_comp_MSE, aes(event,MSE, fill = sampling_design))+
  geom_bar(stat = "identity",position = "dodge")

ggplot(sample_comp_RMSE, aes(event,RMSE, fill = sampling_design))+
  geom_bar(stat = "identity",position = "dodge")+
  ylab("RMSE (m3)")
