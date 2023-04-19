#tracking survival for LC trial
library(janitor)
library(dplyr)

#initial data
initial_trial <- read.csv("2022_growth&inventory_analylsis/growth_raw_collection/LC_growth&inventory_2021/LC_Aug_2021_raw.csv")


initial_count <- table(initial_trial$event)


initial_count <- aggregate(ID ~ event, data = initial_trial, FUN = length)
initial_count <- initial_count[c(3:17), ]

#final data from 2022 season
final_trial <- read.csv("2022_growth&inventory_analylsis/growth_raw_collection/LC_growth&inventory_sheets_2022/LC_9.7.22_DBH_H.csv", header = TRUE)
#final_trial <- final_trial %>% row_to_names(row_number = 1)

final_trial <- subset(final_trial, Height > 0)
final_count <- aggregate(final_trial$ID, by=list(event=final_trial$event), FUN=length)
#final_count <- final_count[c(2:16), ]

survival <- cbind(initial_count, final_count)
total_survival <- (sum(survival$x))/(sum(survival$ID))

survival$percent <- survival$x/survival$ID

write.csv(survival, file = "LC_survival_2022.csv")
