#Cleaning data for Sukyhun's analysis of growth and inventory

library(ggplot2)
library(dplyr)
library(readxl)

#read in previous year's data

growth_dat22 <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")

#Check to see that trees that are supposed to be excluded are not present

excluded_trees_raw <- read_xlsx("LC_2022/2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/LC_all trees to exclude.xlsx",1)
excluded_trees <- as.vector(excluded_trees_raw$`to exclude - dead, replaced, small, diseased`)

growth_dat22 <- filter(growth_dat22, !ID %in% excluded_trees)

#Check to see that CT1 is not present

growth_dat22 <- filter(growth_dat22, !event_short == "CT1")


#Compile 2023 data
LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2"))
growth_dat_may23 <- inner_join(growth_dat22, LC_meta)


growth_dat_may_2023 <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/05_06_DBHH.csv",skip = 5)
growth_dat_may_2023 <- inner_join(growth_dat_may_2023, LC_meta)
growth_dat_may_2023$D664 <- growth_dat_may_2023$DBH
growth_dat_may_2023$H664 <- growth_dat_may_2023$Height*304.8
growth_dat_may_2023 <- subset(growth_dat_may_2023, select = c(row,column,ID,event_short,block, construct,construct2,D664,H664))

growth_dat_aug_2023 <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/8_7_inventory.csv")
growth_dat_aug_2023 <- inner_join(growth_dat_aug_2023, LC_meta)
growth_dat_aug_2023$D757 <- growth_dat_aug_2023$Diameter
growth_dat_aug_2023 <- subset(growth_dat_aug_2023, select = c(row,column,event_short,block,construct,construct2,ID,D757))

growth_dat_sep_2023 <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/9_21_inventory.csv")
growth_dat_sep_2023 <- inner_join(growth_dat_sep_2023, LC_meta)
growth_dat_sep_2023$D801 <- growth_dat_sep_2023$Caliper
growth_dat_sep_2023$DBH801 <- growth_dat_sep_2023$DBH
growth_dat_sep_2023$H801 <- growth_dat_sep_2023$Height*304.8
growth_dat_sep_2023 <- subset(growth_dat_sep_2023, select = c(row,column,event_short,block,construct,construct2,ID,D801,DBH801,H801))

growth_dat23 <- inner_join(growth_dat_may_2023,growth_dat_aug_2023)
growth_dat23 <- inner_join(growth_dat23, growth_dat_sep_2023)

growth_dat <- inner_join(growth_dat22, growth_dat23)

#Exclude border trees

growth_dat <-subset(growth_dat, row != 1)
growth_dat <-subset(growth_dat, column != 30)

####Break down into height and diameter########

growth_dat_ht <- growth_dat[,c(3:18,45,49)]
growth_dat_diam <- growth_dat[,c(3:10,19:26,44,47)]


#########################################Height Growth #################################################

#find trees that show negative ht growth from one measurement period to the next

growth_dat_ht$delta1 <- growth_dat_ht$H144-growth_dat_ht$H49
mean(growth_dat_ht$delta1)

#isolate the trees that have a negative growth.
ht_delta_1 <- subset(growth_dat_ht, growth_dat_ht$delta1 < 0)
delta1 <- c(row.names(ht_delta_1))

#option 2: Or assume no growth when there is negative results. This is to take care of assumed measurement error.
for (i in delta1) {
  print(i)
  growth_dat_ht[delta1, 10] <- growth_dat_ht[delta1, 9]
}

#Then go for remaining intervals

# Days 144 to 299 (nov 2021 to April 2022)
growth_dat_ht$delta2 <- growth_dat_ht$H299-growth_dat_ht$H144
mean(growth_dat_ht$delta2)


ht_delta_2 <- subset(growth_dat_ht, growth_dat_ht$delta2 < 0)
delta2 <- c(row.names(ht_delta_2))

for (i in delta2) {
  print(i)
  growth_dat_ht[delta2, 11] <- growth_dat_ht[delta2, 10]
}

#Days 299 to 335 (may to june)


growth_dat_ht$delta3 <- growth_dat_ht$H335 - growth_dat_ht$H299
mean(growth_dat_ht$delta3)

ht_delta_3 <- subset(growth_dat_ht, growth_dat_ht$delta3 < 0)
delta3 <- c(row.names(ht_delta_3))

for (i in delta3) {
  print(i)
  growth_dat_ht[delta3, 12] <- growth_dat_ht[delta3, 11]
}

#Day 335 to 357 (june to july)

growth_dat_ht$delta4 <- growth_dat_ht$H357 - growth_dat_ht$H335
mean(growth_dat_ht$delta4)

ht_delta_4 <- subset(growth_dat_ht, growth_dat_ht$delta4 < 0)
delta4 <- c(row.names(ht_delta_4))

for (i in delta4) {
  print(i)
  growth_dat_ht[delta4, 13] <- growth_dat_ht[delta4, 12]
} 


#Dat 17 and 18 are candidates to exclude because their estimation in Dec 2021 was so poor

##Day 357 (july) to 385 (august)

growth_dat_ht$delta5 <- growth_dat_ht$H385 - growth_dat_ht$H357
mean(growth_dat_ht$delta5) 

ht_delta_5 <- subset(growth_dat_ht, growth_dat_ht$delta5 < 0)
delta5 <- c(row.names(ht_delta_5))

for (i in delta5) {
  print(i)
  growth_dat_ht[delta5, 14] <- growth_dat_ht[delta5, 13]
} 

### from 384 (august 1 2022) to 419 (september 7 2022)

growth_dat_ht$delta6 <- growth_dat_ht$H421 - growth_dat_ht$H385
mean(growth_dat_ht$delta6)

ht_delta_6 <- subset(growth_dat_ht, growth_dat_ht$delta6 < 0)
delta6 <- c(row.names(ht_delta_6))

for (i in delta6) {
  print(i)
  growth_dat_ht[delta6, 15] <- growth_dat_ht[delta6, 14]
  
}

#from 419 (september 7 22) to 500 (Nov 20 2022)
growth_dat_ht$delta7 <- growth_dat_ht$H497 - growth_dat_ht$H421
mean(growth_dat_ht$delta7)

ht_delta_7 <- subset(growth_dat_ht, growth_dat_ht$delta7 < 0)
delta7 <- c(row.names(ht_delta_7))

for (i in delta7) {
  print(i)
  growth_dat_ht[delta7, 16] <- growth_dat_ht[delta7, 15]
}

#certainly more negative growth seen in small block. This may actually be shrinking due to water loss


#from 497 (Nov 20 2022) to 664 (May 6 2023)
growth_dat_ht$delta8 <- growth_dat_ht$H664 - growth_dat_ht$H497
mean(growth_dat_ht$delta8)

ht_delta_8 <- subset(growth_dat_ht, growth_dat_ht$delta8 < 0 | growth_dat_ht$delta8 > 500)
delta8 <- c(row.names(ht_delta_8))

#fix incorrect reading in Nov 2022 for LCOR-289
growth_dat_ht[114,16] <- growth_dat_ht[114,17]

for (i in delta8) {
  print(i)
  growth_dat_ht[delta8, 17] <- growth_dat_ht[delta8, 16]
}

#from 664 (May 6 2023) to 801 (Sep 20 2023)
growth_dat_ht$delta9 <- growth_dat_ht$H801 - growth_dat_ht$H664
#Need to check out LCOR-469 and LCOR-167
#remove LCOR 068 from analysis. It was topped
growth_dat_ht <- subset(growth_dat_ht, ID != "LCOR-068")
growth_dat_diam <- subset(growth_dat_diam, ID != "LCOR-068")
growth_dat <- subset(growth_dat, ID != "LCOR-068")

mean(growth_dat_ht$delta9)

ht_delta_9 <- subset(growth_dat_ht, growth_dat_ht$delta9 < 0)
delta9 <- c(row.names(ht_delta_9))

for (i in delta9) {
  print(i)
  growth_dat_ht[delta9, 18] <- growth_dat_ht[delta9, 17]
}




###############################Diameter Growth ##############################################################

#find trees that show negative diam growth from one measurement period to the next

growth_dat_diam$delta1 <- growth_dat_diam$D144-growth_dat_diam$D49
mean(growth_dat_diam$delta1)

#isolate the trees that have a negative growth.
diam_delta_1 <- subset(growth_dat_diam, growth_dat_diam$delta1 < 0)
d_delta1 <- c(row.names(diam_delta_1))

#option 2: Or assume no growth when there is negative results. This is to take care of assumed measurement error.
for (i in d_delta1) {
  print(i)
  growth_dat_diam[d_delta1, 10] <- growth_dat_diam[d_delta1, 9]
}

#Then go for remaining intervals

# Days 144 to 298 (nov 2021 to April 2022)
growth_dat_diam$delta2 <- growth_dat_diam$D299-growth_dat_diam$D144
mean(growth_dat_diam$delta2)


diam_delta_2 <- subset(growth_dat_diam, growth_dat_diam$delta2 < 0)
d_delta2 <- c(row.names(diam_delta_2))

for (i in d_delta2) {
  print(i)
  growth_dat_diam[d_delta2, 11] <- growth_dat_diam[d_delta2, 10]
}

#Days 298 to 334 (may to june)


growth_dat_diam$delta3 <- growth_dat_diam$D335 - growth_dat_diam$D299
mean(growth_dat_diam$delta3)

diam_delta_3 <- subset(growth_dat_diam, growth_dat_diam$delta3 < 0)
d_delta3 <- c(row.names(diam_delta_3))

for (i in d_delta3) {
  print(i)
  growth_dat_diam[d_delta3, 12] <- growth_dat_diam[d_delta3, 11]
}

#Day 334 to 356 (june to july)

growth_dat_diam$delta4 <- growth_dat_diam$D357 - growth_dat_diam$D335
mean(growth_dat_diam$delta4)

diam_delta_4 <- subset(growth_dat_diam, growth_dat_diam$delta4 < 0)
d_delta4 <- c(row.names(diam_delta_4))

for (i in d_delta4) {
  print(i)
  growth_dat_diam[d_delta4, 13] <- growth_dat_diam[d_delta4, 12]
} 

##Day 356 (july) to 384 (august)

growth_dat_diam$delta5 <- growth_dat_diam$D385 - growth_dat_diam$D357
mean(growth_dat_diam$delta5) 

diam_delta_5 <- subset(growth_dat_diam, growth_dat_diam$delta5 < 0)
d_delta5 <- c(row.names(diam_delta_5))

for (i in d_delta5) {
  print(i)
  growth_dat_diam[d_delta5, 14] <- growth_dat_diam[d_delta5, 13]
} 

### from 384 (august 1 2022) to 419 (september 7 2022)

growth_dat_diam$delta6 <- growth_dat_diam$D421 - growth_dat_diam$D385
mean(growth_dat_diam$delta6)

diam_delta_6 <- subset(growth_dat_diam, growth_dat_diam$delta6 < 0)
d_delta6 <- c(row.names(diam_delta_6))

for (i in d_delta6) {
  print(i)
  growth_dat_diam[d_delta6, 15] <- growth_dat_diam[d_delta6, 14]
  
}

#from 419 (september 7 22) to 500 (Nov 20 2022)
growth_dat_diam$delta7 <- growth_dat_diam$D497 - growth_dat_diam$D421
mean(growth_dat_diam$delta7)

diam_delta_7 <- subset(growth_dat_diam, growth_dat_diam$delta7 < 0 | growth_dat_diam$delta7 > 10)
d_delta7 <- c(row.names(diam_delta_7))

for (i in d_delta7) {
  print(i)
  growth_dat_diam[d_delta7, 16] <- growth_dat_diam[d_delta7, 15]
}


#from 497 (Nov 20 2022) to 664 (May 6 2023)
growth_dat_diam$delta8 <- growth_dat_diam$D664 - growth_dat_diam$D497
mean(growth_dat_diam$delta8)

d_delta_8 <- subset(growth_dat_diam, growth_dat_diam$delta8 < 0)
d_delta8 <- c(row.names(d_delta_8))

for (i in d_delta8) {
  print(i)
  growth_dat_diam[delta8, 17] <- growth_dat_diam[delta8, 16]
}

#from 664 (May 6 2023) to 801 (Sep 20 2023)
growth_dat_diam$delta9 <- growth_dat_diam$D801 - growth_dat_diam$D664
growth_dat_diam <- na.omit(growth_dat_diam)
mean(growth_dat_diam$delta9)

d_delta_9 <- subset(growth_dat_diam, growth_dat_diam$delta9 < 0)
d_delta9 <- c(row.names(d_delta_9))

for (i in d_delta9) {
  print(i)
  growth_dat_diam[delta9, 18] <- growth_dat_diam[delta9, 17]
}


#There could be a need for height and diameter imputation - either predicted by height to predicted by neighbors.

library(dplyr)
growth_dat_ht_clean <- growth_dat_ht[,-c(19:27)]
growth_dat_diam_clean <- growth_dat_diam[,-c(19:27)]
growth_dat_clean <- inner_join(growth_dat_ht_clean, growth_dat_diam_clean)

##finding outliers

ht_model <- lm(H801~event+block, data = growth_dat_clean)
summary(ht_model)
plot(ht_model, which = 1)
#obs 237, 10 and 17 are candidates
plot(ht_model, which = 2)
plot(ht_model, which = 3)

diam_model <- lm(D801~event+block, data = growth_dat_clean)
summary(diam_model)
plot(diam_model, which = 1)
plot(diam_model, which = 2)
plot(diam_model, which = 3)
plot(diam_model)

#obs 378 is a candidate, standard residual bit greater than 2. Leverage is not high. Will leave in.

ht_model2 <- lm(H49~event+block, data = growth_dat_clean)
summary(ht_model2)
plot(ht_model2, which = 1)
plot(ht_model2, which = 2)
plot(ht_model2, which = 3)
#looks ok

#Calculate volume index

growth_dat_clean$V49 = ((3.14159*(((growth_dat_clean$D49/10)/2)^2))*(growth_dat_clean$H49/10))
growth_dat_clean$V144 = ((3.14159*(((growth_dat_clean$D144/10)/2)^2))*(growth_dat_clean$H144/10))
growth_dat_clean$V299 = ((3.14159*(((growth_dat_clean$D299/10)/2)^2))*(growth_dat_clean$H299/10))
growth_dat_clean$V335 = ((3.14159*(((growth_dat_clean$D335/10)/2)^2))*(growth_dat_clean$H335/10))
growth_dat_clean$V357 = ((3.14159*(((growth_dat_clean$D357/10)/2)^2))*(growth_dat_clean$H357/10))
growth_dat_clean$V385 = ((3.14159*(((growth_dat_clean$D385/10)/2)^2))*(growth_dat_clean$H385/10))
growth_dat_clean$V421 = ((3.14159*(((growth_dat_clean$D421/10)/2)^2))*(growth_dat_clean$H421/10))
growth_dat_clean$V497 = ((3.14159*(((growth_dat_clean$D497/10)/2)^2))*(growth_dat_clean$H497/10))
growth_dat_clean$V664 = ((3.14159*(((growth_dat_clean$D664/10)/2)^2))*(growth_dat_clean$H664/10))
growth_dat_clean$V801 = ((3.14159*(((growth_dat_clean$D801/10)/2)^2))*(growth_dat_clean$H801/10))

write.csv(growth_dat_clean, file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")

