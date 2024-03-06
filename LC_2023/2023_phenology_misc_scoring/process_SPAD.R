#process SPAD data
library(dplyr)
library(tidyverse)


growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")

SPAD_2022_1 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/LC_SPAD_Compilation_07_2022.csv")
SPAD_2022_2 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/LC_August_SPAD.csv")
SPAD_2022_3 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/Sep_SPAD.csv")

SPAD_2023_1 <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/Drought_Scores_SPAD_compiled.csv")

str(SPAD_2022_1)
SPAD_2022_1$SPAD_June_2022 <- as.numeric(SPAD_2022_1$SPAD_June)
SPAD_2022_1$SPAD_July_2022 <- as.numeric(SPAD_2022_1$SPAD_July)
SPAD_2022_1 <- SPAD_2022_1 %>% drop_na(SPAD_June_2022)
SPAD_2022_1 <- SPAD_2022_1 %>% drop_na(SPAD_July_2022)


str(SPAD_2022_2)
SPAD_2022_2$SPAD_aug_2022 <- as.numeric(SPAD_2022_2$Aug_avg_SPAD)
SPAD_2022_2 <- SPAD_2022_2 %>% drop_na(SPAD_aug_2022)

str(SPAD_2022_3)
SPAD_2022_3$SPAD_sep_2022 <- as.numeric(SPAD_2022_3$SPADavg)
SPAD_2022_3 <- SPAD_2022_3 %>% drop_na(SPAD_sep_2022)

str(SPAD_2023_1)

#interpolate missing values with day's average

SPAD_2023_1$SPAD_July_2023 <- as.numeric(SPAD_2023_1$SPAD_18)
SPAD_2023_1_1 <- SPAD_2023_1 %>% drop_na(SPAD_July_2023)
SPAD_july_interpolate_Set <- growth_dat[!growth_dat$ID %in% SPAD_2023_1_1$ID,]
SPAD_2023_1[SPAD_2023_1$ID %in% SPAD_july_interpolate_Set$ID,]$SPAD_July_2023 <- mean(SPAD_2023_1$SPAD_July_2023,na.rm=TRUE)
SPAD_2023_1_1 <- SPAD_2023_1 %>% drop_na(SPAD_July_2023)

SPAD_2023_1$SPAD_Aug_2023 <- as.numeric(SPAD_2023_1$SPAD_32)
SPAD_2023_1_2 <- SPAD_2023_1 %>% drop_na(SPAD_Aug_2023)
SPAD_aug_interpolate_Set <- growth_dat[!growth_dat$ID %in% SPAD_2023_1_2$ID,]
SPAD_2023_1[SPAD_2023_1$ID %in% SPAD_aug_interpolate_Set$ID,]$SPAD_Aug_2023 <- mean(SPAD_2023_1$SPAD_Aug_2023,na.rm=TRUE)


SPAD_2023_1$SPAD_Sep_2023 <- as.numeric(SPAD_2023_1$SPAD_55)
SPAD_2023_1_3 <- SPAD_2023_1 %>% drop_na(SPAD_Sep_2023)
SPAD_sep_interpolate_Set <- growth_dat[!growth_dat$ID %in% SPAD_2023_1_3$ID,]
SPAD_2023_1[SPAD_2023_1$ID %in% SPAD_sep_interpolate_Set$ID,]$SPAD_Sep_2023 <- mean(SPAD_2023_1$SPAD_Sep_2023,na.rm=TRUE)


SPAD_2022_1 <- na.omit(SPAD_2022_1)
SPAD_2022_2 <- na.omit(SPAD_2022_2)
SPAD_2022_3 <- SPAD_2022_3 %>% drop_na(SPAD_sep_2022)


SPAD_2023_1_full <- SPAD_2023_1 %>% drop_na(c(SPAD_July_2023,SPAD_Aug_2023,SPAD_Sep_2023))


SPAD_df_list <- list(SPAD_2022_1,SPAD_2022_2,SPAD_2022_3,SPAD_2023_1)

SPAD_full <- SPAD_df_list %>% reduce(inner_join, by=c("ID"))

SPAD_full[,34]
str(SPAD_full)
SPAD_full$Aug_avg_SPAD <- as.numeric(SPAD_full$Aug_avg_SPAD)
SPAD_summary <- SPAD_full[,c(1,2,3,6,7,15,22,32,33,34)]
str(SPAD_summary)
SPAD_summary$avg_SPAD <- rowMeans(SPAD_summary[,-c(1,2,3)])
SPAD_summary <- na.omit(SPAD_summary)
         
#iput missing/erroneous values
SPAD_summary[which(SPAD_summary$SPAD_June_2022 > 70),]
SPAD_summary[which(SPAD_summary$SPAD_June_2022 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_July_2022 > 70),]
SPAD_summary[which(SPAD_summary$SPAD_July_2022 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_aug_2022 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_aug_2022 > 70),]
SPAD_summary[which(SPAD_summary$SPAD_sep_2022 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_sep_2022 > 70),]
#imput row 372 with mean for that collection
SPAD_summary[372,5] <- mean(SPAD_summary$SPAD_sep_2022)
SPAD_summary[which(SPAD_summary$SPAD_July_2023 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_July_2023 > 70),]
SPAD_summary[which(SPAD_summary$SPAD_Aug_2023 < 10),]
#imput row 7 with mean for that collection
SPAD_summary[7,7] <- mean(SPAD_summary$SPAD_Aug_2023)
SPAD_summary[which(SPAD_summary$SPAD_Aug_2023 > 70),]
SPAD_summary[which(SPAD_summary$SPAD_Sep_2023 < 10),]
SPAD_summary[which(SPAD_summary$SPAD_Sep_2023 > 70),]

SPAD_summary$row <- SPAD_summary$row.x
SPAD_summary$column <- SPAD_summary$column.x
SPAD_summary <- SPAD_summary[,-c(1,2)]
write.csv(SPAD_summary, file = "LC_2023/2023_phenology_misc_scoring/SPAD_summary.csv")

#Explore trends in SPAD data

cormatrix <- (SPAD_summary[,-c(1,9,10,11)])
colnames(cormatrix) <- c("6_22","7_22","8_22","9_22","7_23","8_23","9_23")
#install.packages("corrplot")
library(corrplot)
corrplot.mixed(cor(cormatrix), upper = "shade",lower = "number")
#save
png(filename = "LC_2023/2023_growth_inventory_analysis/figures/SPAD_cor_ful.png")
corrplot.mixed(cor(cormatrix), upper = "shade",lower = "number")
dev.off()

#look at large block only
SPAD_LB <- subset(SPAD_summary,column < 30)
LB_matrix <- SPAD_LB[,-c(1,9,10,11)]
colnames(LB_matrix) <- c("6_22","7_22","8_22","9_22","7_23","8_23","9_23")
corrplot.mixed(cor(LB_matrix), upper = "shade",lower = "number")
#save
png(filename = "LC_2023/2023_growth_inventory_analysis/figures/SPAD_cor_LB.png")
corrplot.mixed(cor(LB_matrix), upper = "shade",lower = "number")
dev.off()
#fairly robust anticorrelation between June and July 2022
#moderate correlation between all months measured in 2023

#look at small block only
SPAD_SB <- subset(SPAD_summary,column > 30)
SB_matrix <- SPAD_SB[,-c(1,9,10,11)]
colnames(SB_matrix) <- c("6_22","7_22","8_22","9_22","7_23","8_23","9_23")
corrplot.mixed(cor(SB_matrix), upper = "shade",lower = "number")

png(filename = "LC_2023/2023_growth_inventory_analysis/figures/SPAD_cor_SB.png")
corrplot.mixed(cor(SB_matrix), upper = "shade",lower = "number")
dev.off()
#Signal:noise stronger
#same anti-correlation seen between June and July 2022,
#correlation between July, August and September seen in both years.


#overall conclusions - June 2022 was weird! not strongly related to anything by anti-correlated with the following month
#Correlations more apparent in the small block but trends are consistent in order between blocks
#correlations between from July through September seen for the most part in both years. Stronger in 2023.
#weak to no correlations from predicting one year to the next.

#from a biological standpoint this may suggest that effects in SPAD are more apparent when there is some environmental heterogeneity
#No pressure for genetic variation to be expressed under healthy, well-watered conditions?

#examine strange relationship between June and July 2022

LC_meta <- read.csv("LC_2023/2023_growth_inventory_analysis/LivingCarbon_9_20_metadata.csv")
SPAD_summary_2 <- inner_join(LC_meta, SPAD_summary, by = "ID")
ggplot(SPAD_summary_2, aes(x=SPAD_June_2022,y=SPAD_July_2022,color=block))+
         geom_point()

#I suppose the interpretation for this is that one can predict SPAD from these variables, just not in the expected way. If did well in June, tended not to score high in July. Curious
  

#best relationship seen
ggplot(subset(SPAD_summary_2,block=="small"), aes(x=SPAD_July_2023,y=SPAD_Sep_2023))+
  geom_point()


#Examine whether SPAD has any potential to explain variation in volume
growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/processed/LC_9_20_growth_data_cleaned.csv")
growth_SPAD <- subset(growth_dat, select = c("ID","V144","V497","V801"))

growth_SPAD_comp <- inner_join(SPAD_summary,growth_SPAD, by = "ID")
growth_SPAD_comp <- growth_SPAD_comp[,-c(10,11)]
str(growth_SPAD_comp)
colnames(growth_SPAD_comp) <- c("ID","6_22","7_22","8_22","9_22","7_23","8_23","9_23","avg","yr1vol","yr2vol","yr3vol")

#simple corr plot
pdf(file = "LC_2023/2023_growth_inventory_analysis/figures/vol_SPAD_correlation.pdf")
corplt_SPAD <- corrplot.mixed(cor(growth_SPAD_comp[,-1]),upper = "shade",lower = "number")
dev.off()


