#Cleaning data for Sukyhun's analysis of growth and inventory

library(ggplot2)
library(dplyr)
library(readxl)

#read in data

growth_dat <- read.csv(file = "2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/3_22_DBH_H_excluded.csv")


#Check to see that trees that are supposed to be excluded are not present

excluded_trees_raw <- read_xlsx("2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/LC_all trees to exclude.xlsx",1)
excluded_trees <- as.vector(excluded_trees_raw$`to exclude - dead, replaced, small, diseased`)

growth_dat <- filter(growth_dat, !ID %in% excluded_trees)

#Check to see that CT1 is not present

growth_dat <- filter(growth_dat, !event_short == "CT1")

####Break down into height and diameter########

growth_dat_ht <- growth_dat[,c(1,2,3,4,5,6,7,8,9,11,15,19,23,27,31,38)]
growth_dat_diam <- growth_dat[,c(1,2,3,4,5,6,7,8,10,12,16,20,24,28,32,37)]


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

diam_delta_7 <- subset(growth_dat_diam, growth_dat_diam$delta7 < 0)
d_delta7 <- c(row.names(diam_delta_7))

for (i in d_delta7) {
  print(i)
  growth_dat_diam[d_delta7, 16] <- growth_dat_diam[d_delta7, 15]
}

#There could be a need for height and diameter imputation - either predicted by height to predicted by neighbors.

library(dplyr)
growth_dat_ht_clean <- growth_dat_ht[,-c(17:23)]
growth_dat_diam_clean <- growth_dat_diam[,-c(17:23)]
growth_dat_clean <- inner_join(growth_dat_ht_clean, growth_dat_diam_clean)

##finding outliers

ht_model <- lm(H497~event+block, data = growth_dat_clean)
summary(ht_model)
plot(ht_model, which = 1)
#obs 378, 10 and 343 are candidates
plot(ht_model, which = 2)
plot(ht_model, which = 3)

diam_model <- lm(D497~event+block, data = growth_dat_clean)
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


write.csv(growth_dat_clean, file = "LC_3_22_growth_data_cleaned.csv")


