#Fv/Fm analysis
library(janitor)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ggsignif)

Fv_fm_dat_8_1 <- read.csv("LC_2023/2023_physiology_analysis/Other_phys/Fv_fm_8_1_23/2023-08-01/Fv_fm_8_1_filled.csv")
Fv_fm_dat_8_2 <- read.csv("LC_2023/2023_physiology_analysis/Other_phys/Fv_fm_8_2_23/2023-08-02/8_2_23_fv_fm_filled.csv")

Fv_fm_dat <-  rbind(Fv_fm_dat_8_1, Fv_fm_dat_8_2)
Fv_fm_dat <- Fv_fm_dat[c(-2),]
Fv_fm_dat <- row_to_names(clean_names(Fv_fm_dat),1)
Fv_fm_dat <- Fv_fm_dat[,c(-7,-8)]

LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2","H497"))


Fv_fm_dat$`Fv/Fm` <- as.numeric(Fv_fm_dat$`Fv/Fm`)

Fv_fm_dat <- Fv_fm_dat %>% group_by("ID")

Fv_fm_summ <- as.data.frame(aggregate(`Fv/Fm` ~ ID, Fv_fm_dat, mean))
Fv_fm_summ <- inner_join(Fv_fm_summ, LC_meta)
Fv_fm_summ <- Fv_fm_summ%>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~ "control",
  event_short == "CT3" ~ "WT"))


Fv_fm_event <- Fv_fm_summ %>% group_by(event_short) %>% dplyr::summarise(
  n = n(),
  Fvfm_sd = sd(`Fv/Fm`),
  Fvfm = mean(`Fv/Fm`),
  Fvfm_se = (Fvfm_sd/(sqrt(n)))
  
)

Fv_fm_tier <- Fv_fm_summ %>% group_by(tier) %>% dplyr::summarise(
  n = n(),
  Fvfm_sd = sd(`Fv/Fm`),
  Fvfm = mean(`Fv/Fm`),
  Fvfm_se = (Fvfm_sd/(sqrt(n)))
  
)

my_colors <- c("grey0","indianred3","chartreuse4","gray")

Fv_fm_event <- Fv_fm_event%>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "elite",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D" ~ "control", 
  event_short == "CT3" ~ "WT"))

ggplot(Fv_fm_event, aes(y = Fvfm, x =event_short, color = tier))+
  geom_point()+
  scale_color_manual(values = my_colors)+
  geom_errorbar(ymin = Fvfm-Fvfm_se, ymax = Fvfm + Fvfm_se)


Fv_fm_plot <- ggplot(Fv_fm_event, aes(x = event_short, y=Fvfm, color = tier))+
  geom_point()+
  geom_errorbar(aes(ymin = Fvfm-2*Fvfm_se, ymax = Fvfm + 2*Fvfm_se))+
  scale_color_manual(values = my_colors)+
  geom_signif(comparisons = list(c("2H","5A")),
              map_signif_level=TRUE, y_position = 0.81)+
  geom_signif(comparisons = list(c("CT3","5A")),
              map_signif_level=TRUE, y_position = 0.815)+
  ylab("Fv/Fm (Photochemical potential)")+
  xlab("Event")+
  theme_bw()
  
  
Fv_fm_plot

Fv_fm_event$Fvfm
ggplot(Fv_fm_tier)

ggsave(filename = "LC_2023/2023_physiology_analysis/Other_phys/July_Fv_fm.png", plot=Fv_fm_plot, dpi =300)


  
geom_errorbar(ymin = Fv_fm-Fv_fm_se, ymax = Fv_fm + Fv_fm_se)
#Evaluate statistical differences


##Read in predawn water potential data
July_WP <- read.csv("LC_2023/2023_physiology_analysis/diurnal_experiments/8_1_23_predawn_WP_filled.csv")
mean(July_WP$Predawn_MPA)
sd(July_WP$Predawn_MPA)/(sqrt(27))

mean(July_WP$Predawn_MPA) - 2*(sd(July_WP$Predawn_MPA)/(sqrt(27)))
mean(July_WP$Predawn_MPA) + 2*(sd(July_WP$Predawn_MPA)/(sqrt(27)))

