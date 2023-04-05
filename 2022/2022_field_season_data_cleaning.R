#Cleaning data for Sukyhun's analysis

#read in data

firstyr_dat <- read.csv(file = "3_22_DBH_H_excluded.csv")


firstyr_dat_ht <- firstyr_dat[,c(1,2,3,4,5,6,7,8,9,11,15,19,23,27,31,38)]
firstyr_dat_diam <- firstyr_dat[,c(1,2,3,4,5,6,7,8,10,12,16,20,24,28,32,37)]


#########################################Height Growth #################################################

#find trees that show negative ht growth from one measurement period to the next

firstyr_dat_ht$delta1 <- firstyr_dat_ht$H144-firstyr_dat_ht$H49
mean(firstyr_dat_ht$delta1)

#isolate the trees that have a negative growth.
ht_delta_1 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta1 < 0)
delta1 <- c(row.names(ht_delta_1))

#option 2: Or assume no growth when there is negative results. This is to take care of assumed measurement error.
for (i in delta1) {
  print(i)
  firstyr_dat_ht[delta1, 10] <- firstyr_dat_ht[delta1, 9]
}

#Then go for remaining intervals

# Days 144 to 299 (nov 2021 to April 2022)
firstyr_dat_ht$delta2 <- firstyr_dat_ht$H299-firstyr_dat_ht$H144
mean(firstyr_dat_ht$delta2)


ht_delta_2 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta2 < 0)
delta2 <- c(row.names(ht_delta_2))

for (i in delta2) {
  print(i)
  firstyr_dat_ht[delta2, 11] <- firstyr_dat_ht[delta2, 10]
}

#Days 299 to 335 (may to june)


firstyr_dat_ht$delta3 <- firstyr_dat_ht$H335 - firstyr_dat_ht$H299
mean(firstyr_dat_ht$delta3)

ht_delta_3 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta3 < 0)
delta3 <- c(row.names(ht_delta_3))

for (i in delta3) {
  print(i)
  firstyr_dat_ht[delta3, 12] <- firstyr_dat_ht[delta3, 11]
}

#Day 335 to 357 (june to july)

firstyr_dat_ht$delta4 <- firstyr_dat_ht$H357 - firstyr_dat_ht$H335
mean(firstyr_dat_ht$delta4)

ht_delta_4 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta4 < 0)
delta4 <- c(row.names(ht_delta_4))

for (i in delta4) {
  print(i)
  firstyr_dat_ht[delta4, 13] <- firstyr_dat_ht[delta4, 12]
} 


#Dat 17 and 18 are candidates to exclude because their estimation in Dec 2021 was so poor

##Day 357 (july) to 385 (august)

firstyr_dat_ht$delta5 <- firstyr_dat_ht$H385 - firstyr_dat_ht$H357
mean(firstyr_dat_ht$delta5) 

ht_delta_5 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta5 < 0)
delta5 <- c(row.names(ht_delta_5))

for (i in delta5) {
  print(i)
  firstyr_dat_ht[delta5, 14] <- firstyr_dat_ht[delta5, 13]
} 

### from 384 (august 1 2022) to 419 (september 7 2022)

firstyr_dat_ht$delta6 <- firstyr_dat_ht$H421 - firstyr_dat_ht$H385
mean(firstyr_dat_ht$delta6)

ht_delta_6 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta6 < 0)
delta6 <- c(row.names(ht_delta_6))

for (i in delta6) {
  print(i)
  firstyr_dat_ht[delta6, 15] <- firstyr_dat_ht[delta6, 14]

}

#from 419 (september 7 22) to 500 (Nov 20 2022)
firstyr_dat_ht$delta7 <- firstyr_dat_ht$H497 - firstyr_dat_ht$H421
mean(firstyr_dat_ht$delta7)

ht_delta_7 <- subset(firstyr_dat_ht, firstyr_dat_ht$delta7 < 0)
delta7 <- c(row.names(ht_delta_7))

for (i in delta7) {
  print(i)
  firstyr_dat_ht[delta7, 16] <- firstyr_dat_ht[delta7, 15]
}

#certainly more negative growth seen in small block. This may actually be shrinking due to water loss


###############################Diameter Growth ##############################################################

#find trees that show negative diam growth from one measurement period to the next

firstyr_dat_diam$delta1 <- firstyr_dat_diam$D144-firstyr_dat_diam$D49
mean(firstyr_dat_diam$delta1)

#isolate the trees that have a negative growth.
diam_delta_1 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta1 < 0)
d_delta1 <- c(row.names(diam_delta_1))

#option 2: Or assume no growth when there is negative results. This is to take care of assumed measurement error.
for (i in d_delta1) {
  print(i)
  firstyr_dat_diam[d_delta1, 10] <- firstyr_dat_diam[d_delta1, 9]
}

#Then go for remaining intervals

# Days 144 to 298 (nov 2021 to April 2022)
firstyr_dat_diam$delta2 <- firstyr_dat_diam$D299-firstyr_dat_diam$D144
mean(firstyr_dat_diam$delta2)


diam_delta_2 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta2 < 0)
d_delta2 <- c(row.names(diam_delta_2))

for (i in d_delta2) {
  print(i)
  firstyr_dat_diam[d_delta2, 11] <- firstyr_dat_diam[d_delta2, 10]
}

#Days 298 to 334 (may to june)


firstyr_dat_diam$delta3 <- firstyr_dat_diam$D335 - firstyr_dat_diam$D299
mean(firstyr_dat_diam$delta3)

diam_delta_3 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta3 < 0)
d_delta3 <- c(row.names(diam_delta_3))

for (i in d_delta3) {
  print(i)
  firstyr_dat_diam[d_delta3, 12] <- firstyr_dat_diam[d_delta3, 11]
}

#Day 334 to 356 (june to july)

firstyr_dat_diam$delta4 <- firstyr_dat_diam$D357 - firstyr_dat_diam$D335
mean(firstyr_dat_diam$delta4)

diam_delta_4 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta4 < 0)
d_delta4 <- c(row.names(diam_delta_4))

for (i in d_delta4) {
  print(i)
  firstyr_dat_diam[d_delta4, 13] <- firstyr_dat_diam[d_delta4, 12]
} 

##Day 356 (july) to 384 (august)

firstyr_dat_diam$delta5 <- firstyr_dat_diam$D385 - firstyr_dat_diam$D357
mean(firstyr_dat_diam$delta5) 

diam_delta_5 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta5 < 0)
d_delta5 <- c(row.names(diam_delta_5))

for (i in d_delta5) {
  print(i)
  firstyr_dat_diam[d_delta5, 14] <- firstyr_dat_diam[d_delta5, 13]
} 

### from 384 (august 1 2022) to 419 (september 7 2022)

firstyr_dat_diam$delta6 <- firstyr_dat_diam$D421 - firstyr_dat_diam$D385
mean(firstyr_dat_diam$delta6)

diam_delta_6 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta6 < 0)
d_delta6 <- c(row.names(diam_delta_6))

for (i in d_delta6) {
  print(i)
  firstyr_dat_diam[d_delta6, 15] <- firstyr_dat_diam[d_delta6, 14]
  
}

#from 419 (september 7 22) to 500 (Nov 20 2022)
firstyr_dat_diam$delta7 <- firstyr_dat_diam$D497 - firstyr_dat_diam$D421
mean(firstyr_dat_diam$delta7)

diam_delta_7 <- subset(firstyr_dat_diam, firstyr_dat_diam$delta7 < 0)
d_delta7 <- c(row.names(diam_delta_7))

for (i in d_delta7) {
  print(i)
  firstyr_dat_diam[d_delta7, 16] <- firstyr_dat_diam[d_delta7, 15]
}

#There could be a need for height and diameter imputation - either predicted by height to predicted by neighbors.

library(dplyr)
firstyr_dat_ht_clean <- firstyr_dat_ht[,-c(17:23)]
firstyr_dat_diam_clean <- firstyr_dat_diam[,-c(17:23)]
firstyr_dat_clean <- inner_join(firstyr_dat_ht_clean, firstyr_dat_diam_clean)

##finding outliers

ht_model <- lm(H497~event+block, data = firstyr_dat_clean)
summary(ht_model)
plot(ht_model, which = 1)
#obs 378, 10 and 343 are candidates
plot(ht_model, which = 2)
plot(ht_model, which = 3)

diam_model <- lm(D497~event+block, data = firstyr_dat_clean)
summary(diam_model)
plot(diam_model, which = 1)
plot(diam_model, which = 2)
plot(diam_model, which = 3)
plot(diam_model)

#obs 378 is a candidate, standard residual bit greater than 2. Leverage is not high. Will leave in.

ht_model2 <- lm(H49~event+block, data = firstyr_dat_clean)
summary(ht_model2)
plot(ht_model2, which = 1)
plot(ht_model2, which = 2)
plot(ht_model2, which = 3)
#looks ok


write.csv(firstyr_dat_clean, file = "LC_3_22_growth_data_cleaned.csv")


