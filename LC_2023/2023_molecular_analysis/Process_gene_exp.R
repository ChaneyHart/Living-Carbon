#2023 gene expression processing

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
#install.packages("PCAtools")
#library(PCAtools)
library(GGally)
library(tidyverse)

#read in 2021 and 2022 data
Gene_exp_2021_22 <- read.csv("LC_2022/2022_molecular_analysis/2022_gene_expression/Gene_exp_bytree_both_yrs.csv")
#split into 2021 an 22
Gene_exp_21 <- (Gene_exp_2021_22[,c(1:8)])
Gene_exp_21$PLGG1 <- as.numeric(Gene_exp_21$PLGG1_2021)
Gene_exp_21$GDH <- as.numeric(Gene_exp_21$GDH_2021)
Gene_exp_21$MS <- as.numeric(Gene_exp_21$MS_2021)
Gene_exp_21 <- na.omit(Gene_exp_21)
Gene_exp_21$Event <- Gene_exp_21$Event_short
Gene_exp_21 <- Gene_exp_21[,c(1,2,9,10,11)]
Gene_exp_21$sample_date <- "09-12-2021"
Gene_exp_21 <- gather(Gene_exp_21, key = "gene",value = "expression_level",PLGG1,GDH,MS)


Gene_exp_22 <- (Gene_exp_2021_22[,c(1:5,9,10,11)])
Gene_exp_22$PLGG1 <- as.numeric(Gene_exp_22$PLGG1_2022)
Gene_exp_22$GDH <- as.numeric(Gene_exp_22$GDH_2022)
Gene_exp_22$MS <- as.numeric(Gene_exp_22$MS_2022)
Gene_exp_22 <- na.omit(Gene_exp_22)
Gene_exp_22$Event <- Gene_exp_22$Event_short
Gene_exp_22 <- Gene_exp_22[,c(1,2,9,10,11)]
Gene_exp_22$sample_date <- "06-20-2022"
Gene_exp_22 <- gather(Gene_exp_22, key = "gene",value = "expression_level",PLGG1,GDH,MS)

#read in 2023 data
#July
PLGG1_july_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/July_2023_PLGG1.csv")
PLGG1_july_23$ID <- substr(PLGG1_july_23$ID,1,8)
PLGG1_july_23$PLGG1<- PLGG1_july_23$Fold
PLGG1_july_23 <- subset(PLGG1_july_23, select = c(ID,Event,PLGG1))


GDH_july_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/July_2023_GDH.csv")
GDH_july_23$ID <- substr(GDH_july_23$ID,1,8)
GDH_july_23$GDH <- GDH_july_23$Fold
GDH_july_23 <- subset(GDH_july_23, select = c(ID,Event,GDH))


MS_july_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/July_molecular/July_2023_MS.csv")
MS_july_23$ID <- substr(MS_july_23$ID,1,8)
MS_july_23$MS <- MS_july_23$Fold
MS_july_23 <- subset(MS_july_23, select = c(ID,Event,MS))


Gene_exp_July_23 <- inner_join(PLGG1_july_23,GDH_july_23)
Gene_exp_July_23 <- inner_join(Gene_exp_July_23,MS_july_23)
Gene_exp_July_23$sample_date <- "07-14-2023"

Gene_exp_July_23 <- gather(Gene_exp_July_23, key = "gene",value = "expression_level",PLGG1,GDH,MS)
Gene_exp_July_23[142,5] <- 0

#August

PLGG1_aug_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/August_molecular/Aug_2023_PLGG1.csv")
PLGG1_aug_23$ID <- substr(PLGG1_aug_23$ID,1,8)
PLGG1_aug_23$PLGG1 <- PLGG1_aug_23$Fold
PLGG1_aug_23 <- subset(PLGG1_aug_23, select = c(ID,Event,PLGG1))


GDH_aug_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/August_molecular/Aug_2023_GDH.csv")
GDH_aug_23$ID <- substr(GDH_aug_23$ID,1,8)
GDH_aug_23$GDH <- GDH_aug_23$Fold
GDH_aug_23 <- subset(GDH_aug_23, select = c(ID,Event,GDH))


MS_aug_23 <- read.csv(file = "LC_2023/2023_molecular_analysis/August_molecular/Aug_2023_MS.csv")
MS_aug_23$ID <- substr(MS_aug_23$ID,1,8)
MS_aug_23$MS <- MS_aug_23$Fold
MS_aug_23 <- subset(MS_aug_23, select = c(ID,Event,MS))

Gene_exp_aug_23 <- inner_join(PLGG1_aug_23,GDH_aug_23)
Gene_exp_aug_23 <- inner_join(Gene_exp_aug_23,MS_aug_23)
Gene_exp_aug_23$sample_date <- "08-14-2023"

Gene_exp_aug_23 <- gather(Gene_exp_aug_23, key = "gene",value = "expression_level",PLGG1,GDH,MS)
#imput 0s for missing values for controls
Gene_exp_aug_23[133,5] <- 0
Gene_exp_aug_23[135,5] <- 0
Gene_exp_aug_23[136,5] <- 0

#Combine all datasets
str(Gene_exp_21)
str(Gene_exp_22)
str(Gene_exp_July_23)
Gene_exp_July_23$expression_level <- as.numeric(Gene_exp_July_23$expression_level)
Gene_exp_aug_23$expression_level <- as.numeric(Gene_exp_aug_23$expression_level)

gene_exp <- rbind(Gene_exp_21,Gene_exp_22,Gene_exp_July_23,Gene_exp_aug_23)

#imput missing values

missing_names <- subset(gene_exp, Event == "#N/A")
missing_names_2 <- subset(missing_names, ID == "LCOR-057")
replacements <- c(rownames(missing_names_2))

for (i in replacements) {
  print(i)
  gene_exp[replacements,2 ] <- "2H"
}

duplicates <- subset(missing_names, ID != "LCOR-057")
duplicates_list <- c(rownames(duplicates))
#omit duplicates

gene_exp <- gene_exp[-c(370,495,620),]

gene_exp <- gene_exp %>% mutate(Class = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "intermediate",
  Event == "13-15E" | Event == "2H" ~ "high",
  Event == "16-20" | Event== "8-9D" ~ "Control",
  Event == "CT3" ~ "WT"))

gene_exp_full <- na.omit(gene_exp)
write.csv(file = "LC_2023/2023_molecular_analysis/gene_exp.csv",gene_exp_full)

gene_exp_summary <- gene_exp_full %>% group_by(Event,gene) %>% summarize(
  n = n(),
  expression_level_sd = sd(expression_level),
  expression_level_mean = mean(expression_level),
  expression_level_se = expression_level_sd/(sqrt(n))
  
)

gene_exp_summary <- gene_exp_summary %>% mutate(Class = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "intermediate",
  Event == "13-15E" | Event == "2H" ~ "high",
  Event == "16-20" | Event== "8-9D" ~ "Control",
  Event == "CT3" ~ "WT"))

my_colors <- c("gray0","chartreuse4","indianred3","grey")

gene_exp_summary_plot <- ggplot(subset(gene_exp_summary,Class != "WT"),aes(x=Event,y=expression_level_mean,color=Class))+
  geom_point()+
  geom_errorbar(aes(ymax=expression_level_mean + 2*(expression_level_se),ymin=expression_level_mean - 2*(expression_level_se)))+
  facet_wrap( ~ gene)+
  scale_color_manual(values = my_colors)+
  ylab("Relative gene expression level")+
  theme_bw()

gene_exp_summary_plot
ggsave(file="LC_2023/2023_molecular_analysis/gene_exp_summary_plot.png",plot=gene_exp_summary_plot,height = 6, width = 10, units = "in",dpi = 300)  

gene_exp_summary_yearly <- gene_exp_full %>% group_by(Event,gene,sample_date) %>% summarize(
  n = n(),
  expression_level_sd = sd(expression_level),
  expression_level_mean = mean(expression_level),
  expression_level_se = expression_level_sd/(sqrt(n))
  
)

gene_exp_summary_yearly <- gene_exp_summary_yearly %>% mutate(Class = case_when(
  Event == "5A" | Event == "5C" | Event == "4A" ~ "intermediate",
  Event == "13-15E" | Event == "2H" ~ "high",
  Event == "16-20" | Event== "8-9D" ~ "Control",
  Event == "CT3" ~ "WT"))

ggplot(gene_exp_summary_yearly,aes(x=Event,y=expression_level_mean,color=gene))+
  geom_point()+
  geom_errorbar(aes(ymax=expression_level_mean + 2*(expression_level_se),ymin=expression_level_mean - 2*(expression_level_se)))+
  facet_wrap( ~ sample_date)


weather <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")



weather$Datetime <- ymd_hms(weather$Datetime)
weather$Date <- as.Date(weather$Datetime)
weather$Date
weather_8_14 <- subset(weather, Date == "2023-08-14")
weather_8_14$hour <- hour(weather_8_14$Datetime)

weather_long <- gather(weather, key = "meas",value="value",air_temp_max,VPD_max,solar_radiation_max)

ggplot(subset(weather_long, Date > "06-01-2024" & Date < "09-01-2024"),aes(x=Date))+
  geom_point(aes(y=value))+
  geom_line(aes(y=value))+
  facet_wrap(~ meas)

