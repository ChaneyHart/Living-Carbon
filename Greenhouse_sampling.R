#load packages needed
library(dplyr)

#import inventory file
greenhouse_inv <- read.csv(file="LC_2023/Greenhouse_study/Greenhouse_inventory.csv")

greenhouse_inv <- subset(greenhouse_inv, Class != "WT")

set.seed(45)
greenhouse_sample <- greenhouse_inv %>% group_by(Class) %>% sample_n(4)
str(greenhouse_sample)

greenhouse_sample[,]
set.seed(48)
greenhouse_sample_shuffled <- greenhouse_sample[sample(1:nrow(greenhouse_sample)),]

write.csv(greenhouse_sample, file = "LC_2023/Greenhouse_study/greenhouse_phys_sample.csv")
