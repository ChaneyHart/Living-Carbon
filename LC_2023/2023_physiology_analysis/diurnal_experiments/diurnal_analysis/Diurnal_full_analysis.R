##Combined analysis of diurnal experiments

library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(tidyr)
library(stringi)
library(janitor)

#June

June_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_tree.csv)")
June_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_tree.csv)")

June_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_event.csv)")
June_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_event.csv)")

June_summary_event <- inner_join(June_summary_6800_event, June_summary_600_event, by = c("event_short","Timepoint"))

June_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_6800_diurnal_summary_class.csv)")
June_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_26_600_diurnal_summary_class.csv)")

June_summary_class <- inner_join(June_summary_6800_class, June_summary_600_class, by = c("Class","Timepoint"))


#July

July_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_tree.csv)")
July_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_tree.csv)")

July_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_event.csv)")
July_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_event.csv)")

July_summary_event <- inner_join(July_summary_6800_event, July_summary_600_event, by = c("event_short","Timepoint","treatment"))

July_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_6800_diurnal_summary_class.csv)")
July_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_600_diurnal_summary_class.csv)")

July_summary_class <- inner_join(July_summary_6800_class, July_summary_600_class, by = c("Class","Timepoint","treatment"))
July_fluor_set_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_event.csv)")
July_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/July_26_fluor_set_diurnal_summary_class.csv)")

#August

August_summary_6800_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_tree.csv)")
August_summary_600_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_tree.csv)")
August_predawn_summary_tree <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_tree.csv)")

August_summary_6800_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_event.csv)")
August_summary_600_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_event.csv)")
August_predawn_summary_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_event.csv)")

August_summary_event <- inner_join(August_summary_6800_event, August_summary_600_event, by = c("event_short","Timepoint","treatment"))

August_summary_6800_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_6800_diurnal_summary_class.csv)")
August_summary_600_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_summary_class.csv)")
August_predawn_summary_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_600_diurnal_predawn_summary_class.csv)")

August_summary_class <- inner_join(August_summary_6800_class, August_summary_600_class, by = c("Class","Timepoint","treatment"))

August_fluor_set_event <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_event.csv)")
August_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/August_28_fluor_set_diurnal_summary_class.csv)")
#Plot

#define colors
library(RColorBrewer)
my_colors <- c("gray0","indianred3","chartreuse4")
my_colors2 <- c("chartreuse4","gray0","chartreuse4","indianred3","indianred3","indianred3","gray0","gray0")

library(ggsignif)
library(gridExtra)
July_assimilation <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)

July_assimilation

July_gsw <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")
July_gsw

July_summary_class$WUE <- July_summary_class$A/July_summary_class$gsw

July_WUE <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = WUE), size = 2)+
  geom_line(aes(y = WUE, linetype = treatment))+
  scale_color_manual(values = my_colors)+
  ylab("Water use efficiency")
July_WUE


July_ETR <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")
July_ETR

July_diurnal_phys_plot <- grid.arrange(July_assimilation, July_gsw, July_ETR, nrow=3)

July_Jo <- ggplot(subset(July_fluor_set_class), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo), size = 2)+
  geom_line(aes(y = Jo, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo-Jo_se, ymax = Jo+Jo_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")
July_Jo

July_Jo_ratio <- ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  geom_line(aes(y = Jo_ratio, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio+Jo_ratio_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
July_Jo_ratio

#gsw curves

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation")

ggplot(subset(July_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("ETR")




#August
str(August_summary_class)
August_summary_class$Timepoint <- as.numeric(August_summary_class$Timepoint)

August_assimilation_class <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A,linetype = treatment,linetype = treatment))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)

August_assimilation_class

August_assimilation_event <- ggplot(August_summary_event, aes(x = Timepoint, color = event_short, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A,linetype = treatment,linetype = treatment))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors2)+
  ylab("Assimilation (µmol/m^2*s)")

August_assimilation_event

August_gsw <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw,linetype = treatment))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")
August_gsw

August_ETR <- ggplot(August_summary_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")
August_ETR

August_diurnal_phys_plot <- grid.arrange(August_assimilation, August_gsw, August_ETR, nrow=3)

August_Jo <- ggplot(August_fluor_set_class, aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo), size = 2)+
  geom_line(aes(y = Jo,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo-Jo_se, ymax = Jo+Jo_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")
August_Jo

August_Jo_ratio <- ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  geom_line(aes(y = Jo_ratio,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio-Jo_ratio_se, ymax = Jo_ratio+Jo_ratio_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")
August_Jo_ratio

#gsw curves

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = A), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation")

ggplot(subset(August_fluor_set_class, Timepoint != 18), aes(x = gsw, color = Class, shape = treatment))+
  geom_point(aes(y = ETR), size = 2)+
  scale_color_manual(values = my_colors)+
  ylab("ETR")


ggplot(subset(August_fluor_set_class, Class = Elite), aes(x = Timepoint, color = treatment))+
  geom_point(aes(y=A))+
  geom_line(aes(y=A))+
  geom_point(aes(y=ETR))+
  geom_line(aes(y=ETR))

ggplot(subset(August_fluor_set_class, Class = Elite), aes(x = Timepoint, color = treatment))+
  geom_point(aes(y=Jo))+
  geom_line(aes(y=Jo))+
  geom_point(aes(y=Jc))+
  geom_line(aes(y=Jc))
