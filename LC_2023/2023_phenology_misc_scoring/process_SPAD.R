#process SPAD data

growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")

SPAD_2022_1 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/LC_SPAD_Compilation_07_2022.csv")
SPAD_2022_2 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/LC_August_SPAD.csv")
SPAD_2022_3 <- read.csv(file = "LC_2022/2022_leaf_trait_analysis/SPAD/Sep_SPAD.csv")

SPAD_2023_1 <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/Drought_Scores_SPAD_compiled.csv")

str(SPAD_2022_1)
SPAD_2022_1$SPAD_June_2022 <- as.numeric(SPAD_2022_1$SPAD_June)
SPAD_2022_1$SPAD_July_2022 <- as.numeric(SPAD_2022_1$SPAD_July)

str(SPAD_2022_2)
SPAD_2022_2$SPAD_aug_2022 <- as.numeric(SPAD_2022_2$Aug_avg_SPAD)

str(SPAD_2022_3)
SPAD_2022_3$SPAD_sep_2022 <- as.numeric(SPAD_2022_3$SPADavg)

str(SPAD_2023_1)

SPAD_2023_1$SPAD_July_2023 <- as.numeric(SPAD_2023_1$SPAD_18)
SPAD_2023_1$SPAD_Aug_2023 <- as.numeric(SPAD_2023_1$SPAD_32)
SPAD_2023_1$SPAD_Sep_2023 <- as.numeric(SPAD_2023_1$SPAD_55)

SPAD_2022_1 <- na.omit(SPAD_2022_1)
SPAD_2022_2 <- na.omit(SPAD_2022_2)
SPAD_2022_3 <- SPAD_2022_3 %>% drop_na(SPAD_sep_2022)
df %>% drop_na(col1)

SPAD_2023_1 <- SPAD_2023_1 %>% drop_na(c(SPAD_July_2023,SPAD_Aug_2023,SPAD_Sep_2023))


SPAD_df_list <- list(SPAD_2022_1,SPAD_2022_2,SPAD_2022_3,SPAD_2023_1)

SPAD_full <- SPAD_df_list %>% reduce(inner_join, by=c("ID"))

SPAD_full[,34]
str(SPAD_full)
SPAD_full$Aug_avg_SPAD <- as.numeric(SPAD_full$Aug_avg_SPAD)
SPAD_full$avg_SPAD <- rowMeans(SPAD_full[,c(6,7,13,15,22,32,33,34)])
         
write.csv(SPAD_full, file = "LC_2023/2023_phenology_misc_scoring/SPAD_full.csv")
