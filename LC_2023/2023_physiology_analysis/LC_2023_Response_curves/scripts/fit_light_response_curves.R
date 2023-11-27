#process fits for A-ci curves and output them in a csv

library(plantecophys)
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(plantecophys)
library(mgcv)

#read in compiled datasets

AQ_compiled <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/light_response_curves_compiled.csv")

#find outliers
outlier_check <- lm(A~Qin,AQ_compiled)
plot(outlier_check, which = 1)
plot(outlier_check, which = 2)
plot(outlier_check, which = 3)

LCOR307 = subset(AQ_compiled, ID == "LCOR-307")
LCOR307_1 = LCOR307[1:8,]
LCOR307_2 = LCOR307[9:15,]
ggplot(LCOR307_1, aes(x=Qin,y=A))+
  geom_point()


LC_meta <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")

LC_meta <- subset(LC_meta, select = c("ID","event_short","construct","construct2","block"))
AQ_compiled <- inner_join(AQ_compiled, LC_meta)
AQ_compiled <- AQ_compiled %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D"~"escape",
  event_short == "CT3"~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))

AQ_compiled$Q_point <- round(AQ_compiled$Qin, digits=-2)

AQ_tier_summary <- AQ_compiled %>% group_by(tier, Q_point) %>% summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  qt = qt(0.975,n),
  upper_A = A+(qt*A_se),
  lower_A = A-(qt*A_se),
  Qin = mean(Qin)
  
)

AQ_day_summary <- AQ_compiled %>% group_by(sample_date, Q_point) %>% summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  qt = qt(0.975,n),
  upper_A = A+(qt*A_se),
  lower_A = A-(qt*A_se),
  Qin = mean(Qin)
  
)

ggplot(AQ_day_summary, aes(x=Qin,y=A, color = sample_date))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_A,ymax=upper_A))

ggplot(AQ_tier_summary, aes(x=Qin,y=A, color = tier))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_A,ymax=upper_A))
  
