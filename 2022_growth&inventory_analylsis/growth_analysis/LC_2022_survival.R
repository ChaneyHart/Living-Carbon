#tracking survival for LC trial
library(janitor)

#initial data
initial_trial <- read.csv("2022_growth&inventory_analylsis/growth_raw_collection/LC_growth&inventory_2021/LC_Aug_2021_raw.csv")


initial_count <- count(initial_trial$event)
initial_count <- subset(initial_count, x == c("CT717 1", "CT717 3", "LC-102 1", "LC-102 13-15B", "LC-102 13-15E", "LC-102 16-20", "LC-102 1C", "LC-102 2H", "LC-102 4A", "LC-102 4B", "LC-102 5", "LC-102 5A", "LC-102 5C", "LC-102 7", "LC-102 8-9D"))

#final data from 2022 season
final_trial <- read.csv("2022_growth&inventory_analylsis/growth_raw_collection/LC_growth&inventory_sheets_2022/LC_11.20.22_DBH_H.csv", header = TRUE)
final_trial <- final_trial %>% row_to_names(row_number = 1)

final_count <- count(final_trial$event)

survival <- cbind(initial_count, final_count)
