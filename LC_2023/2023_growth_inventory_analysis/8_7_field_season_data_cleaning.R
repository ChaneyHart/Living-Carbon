#Cleaning data for analysis of growth and inventory

library(ggplot2)
library(dplyr)
library(readxl)

#read in data

growth_dat <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/8_7_inventory.csv")


#Check to see that trees that are supposed to be excluded are not present

excluded_trees_raw <- read_xlsx("LC_2022/2022_growth&inventory_analylsis/growth_raw_collection/2022_compiled_growth&inventory/LC_all trees to exclude.xlsx",1)
excluded_trees <- as.vector(excluded_trees_raw$`to exclude - dead, replaced, small, diseased`)

growth_dat <- filter(growth_dat, !ID %in% excluded_trees)

growth_dat <- filter(growth_dat, Diameter > 0)

growth_dat <- filter(growth_dat, Diameter != "bag")

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))

LC_meta <- LC_meta %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

growth_dat <- inner_join(LC_meta, growth_dat, by = "ID")

#old_growth <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
#growth_dat <- inner_join(growth_dat, old_growth, by = "ID")

#Check to see that CT1 is not present

growth_dat <- filter(growth_dat, !event_short == "CT1")
str(growth_dat)
growth_dat$Diameter <- as.numeric(growth_dat$Diameter)
growth_dat$SPAD <- as.numeric(growth_dat$SPAD)
growth_dat$Drought.Score <- as.factor(growth_dat$Drought.Score)

##graph diameter

July_2023_diam_comp <- ggplot(growth_dat, aes(x = reorder(event_short,Diameter), y = Diameter, fill = construct))+
  geom_boxplot()+
  xlab("Event")+
  ylab("July 2023 diameter(mm)")

July_2023_diam_comp

July_2023_summ_event <- growth_dat %>% group_by(event_short) %>% dplyr::summarise(
  n = n(),
  D_sd = sd(Diameter),
  D = mean(Diameter),
  D_se = (D_sd/(sqrt(n))),
)



July_2023_summ_event <- July_2023_summ_event %>% mutate(Class = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" | event_short == "CT3" ~ "control",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "unanalyzed"))

July_2023_summ_event <- July_2023_summ_event %>% mutate(construct = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" | event_short == "13-15E" | event_short == "2H"| event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B" ~ "transgenic",
  event_short == "CT3" ~ "control",
  event_short == "16-20" | event_short == "8-9D" ~ "escape"))

  


my_colors2 <- c("gray0","indianred3","chartreuse4","grey")

July_2023_diam_plot <- ggplot(July_2023_summ_event, aes(x=reorder(event_short, D), y = D, shape = construct, color = Class))+
  geom_point(size=4)+
  xlab("Event")+
  ylab("July 2023 diameter(mm)")+
  geom_errorbar(aes(ymin = D-D_se, ymax = D+D_se), width = 1)+
  scale_color_manual(values = my_colors2)
  

July_2023_diam_plot
ggsave(filename = "July_2023_diam_plot.png", plot = July_2023_diam_plot, height = 8, width = 12, units = "in", dpi=300)

# Additional libraries
library(emmeans)
library(nlme)
library(lme4)

unique(growth_dat$event.x)

# Combine 16-20 and 8-9D
growth_dat$event2 <- growth_dat$event.x
growth_dat$event2[growth_dat$event.x == "LC-102 16-20" |
                growth_dat$event.x == "LC-102 8-9D"] <- "escape"


mod00 <- lm(Diameter ~ D49 + construct2.x + block.x, data = growth_dat)
summary(mod00)

# Make a block as a random effect
mod2 <- lme(Diameter ~ D49 + construct2.x, random = ~1|block.x, data = growth_dat)
summary(mod2)

# Residual and qq plots
plot(fitted(mod2), residuals(mod2), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(mod2)); qqline(residuals(mod2))

#event effect
d_mod4 <- lme(Diameter ~ D49 + event2, random = ~1|block.x, data = growth_dat)
summary(d_mod4)

# Residual and qq plots
plot(fitted(d_mod4), residuals(d_mod4), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(d_mod4)); qqline(residuals(d_mod4))

emmeans(d_mod4, specs = pairwise ~event2)

#Effect size
predict_set_d_II <- subset(growth_dat, select = c("ID","event.x","event_short.x","construct.x","construct2.x","block.x","D49","Diameter"))
predict_set_d_II$predicted_D <- predict(d_mod4, new_data = predict_set_d_II, interval = "confidence")

#average modeled height
modelled_CT3_mean_d <- mean(predict_set_d_II$predicted_D)
modelled_pooled_escape_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.3197
modelled_5A_mean_d <- mean(predict_set_d_II$predicted_D497) + 3.0533
modelled_5C_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.9791
modelled_5_mean_d <- mean(predict_set_d_II$predicted_D497) + 1.5744
modelled_2H_mean_d <- mean(predict_set_d_II$predicted_D497) - 0.3084
modelled_13_15E_mean_d <- mean(predict_set_d_II$predicted_D497) + 0.4132
modelled_13_15B_mean_d <- mean(predict_set_d_II$predicted_D497) + 0.5951

event_diameter_effect <- as.data.frame(c(modelled_CT3_mean_d, modelled_pooled_escape_mean_d, modelled_5A_mean_d, modelled_5C_mean_d, modelled_5_mean_d, modelled_2H_mean_d, modelled_13_15E_mean_d, modelled_13_15B_mean_d))

row.names(event_diameter_effect) <- c("control","escape","5A","5C","5","2H","13_15E","13_15B")
colnames(event_diameter_effect) <- c("diameter")

growth_dat$Drought.Score <- as.numeric(growth_dat$Drought.Score)

#need to convert to wide format for drought scores
library(dplyr)
library(tidyr)
??pivot_wider
growth_dat
drought_dat <- subset(growth_dat, select = c(ID,Drought.Score))
growth_dat$d0 <- ifelse(growth_dat$Drought.Score == "0", 1, 0)
growth_dat$d1 <- ifelse(growth_dat$Drought.Score == "1", 1, 0)
growth_dat$d2 <- ifelse(growth_dat$Drought.Score == "2", 1, 0)
growth_dat$d3 <- ifelse(growth_dat$Drought.Score == "3", 1, 0)
growth_dat$d4 <- ifelse(growth_dat$Drought.Score == "4", 1, 0)
growth_dat$d5 <- ifelse(growth_dat$Drought.Score == "5", 1, 0)

growth_dat_drought_summary_8_7 <- growth_dat %>% group_by(event_short) %>% summarise(
  d0 = sum(d0),
  d1 = sum(d1),
  d2 = sum(d2),
  d3 = sum(d3),
  d4 = sum(d4),
  d5 = sum(d4)
)

ggplot(growth_dat, aes(x = event_short, y = Drought.Score))+
  geom_boxplot()




#graph drought scores
install.packages("tidymodels")
growth_dat_recipe
