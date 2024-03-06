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
June_fluor_set_class <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_analysis/June_28_fluor_set_diurnal_summary_class.csv)")
June_fluor_set_clean <- read.csv(file = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_fluor_set_clean.csv")


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


#######graph######


#define colors
library(RColorBrewer)
my_colors <- c("gray0","chartreuse4","indianred3","grey")

library(ggsignif)
library(gridExtra)
June_assimilation <- ggplot(subset(June_summary_class, Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Assimilation (µmol/m^2*s)")+
  xlab("Hour")+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))


June_assimilation

June_gsw <- ggplot(subset(June_summary_class, Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

June_gsw

June_ETR <- ggplot(subset(June_summary_class, Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR), size = 2)+
  geom_line(aes(y = ETR))+ 
  geom_errorbar(aes(ymin = ETR-ETR_se, ymax = ETR+ETR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))
    
June_ETR

June_ETR_6800 <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR_6800), size = 2)+
  geom_line(aes(y = ETR_6800))+ 
  geom_errorbar(aes(ymin = ETR_6800-ETR_6800_se, ymax = ETR_6800+ETR_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))
    
June_ETR_6800

June_PhiPS2 <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PhiPS2), size = 2)+
  geom_line(aes(y = PhiPS2))+ 
  geom_errorbar(aes(ymin = PhiPS2-PhiPS2_se, ymax = PhiPS2+PhiPS2_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("efficiency of PSII")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

June_PhiPS2

June_PhiPS2_6800 <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PhiPS2_6800), size = 2)+
  geom_line(aes(y = PhiPS2_6800))+ 
  geom_errorbar(aes(ymin = PhiPS2_6800-PhiPS2_6800_se, ymax = PhiPS2_6800+PhiPS2_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("efficiency of PSII")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

June_PhiPS2_6800

June_Tleaf <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = Tleaf), size = 2)+
  geom_line(aes(y = Tleaf))+ 
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

June_Tleaf

June_PAR <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PAR), size = 2)+
  geom_line(aes(y = PAR))+ 
  geom_errorbar(aes(ymin = PAR-PAR_se, ymax = PAR+PAR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*sec)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

June_PAR

June_Ci <- ggplot(subset(June_summary_class,Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = Ci), size = 2)+
  geom_line(aes(y = Ci))+ 
  geom_errorbar(aes(ymin = Ci-Ci_se, ymax = Ci+Ci_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*sec)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

June_Ci

June_diurnal_phys_plot <- grid.arrange(June_assimilation, June_gsw,June_Tleaf,nrow=1)
June_diurnal_light_plot <- grid.arrange(June_PAR, June_PhiPS2_6800, June_ETR_6800,nrow=1)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_diurnal_phys_plot.png", plot = June_diurnal_phys_plot, height = 3, width = 7, units = "in", dpi = 300)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_diurnal_light_plot.png", plot = June_diurnal_light_plot, height = 3, width = 7, units = "in", dpi = 300)


#photorespiration phys

#png(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/A_ETR.png",height = 4, width = 4, units = "in",res = 300)
ggplot(subset(June_fluor_set,Class != "WT"), aes(x=A-R,y=ETR,color = Class))+
  geom_point()+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Gross assimilation (µmol/m^2*s)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

#dev.off()


June_Jo <- ggplot(subset(June_fluor_set_class,Class != "WT"), aes(x = Timepoint, group=Class,color = Class))+
  geom_point(aes(y = Jo_6800), size = 2)+
  geom_line(aes(y = Jo_6800))+
  geom_errorbar(aes(ymin = Jo_6800-Jo_6800_se, ymax = Jo_6800+Jo_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))
June_Jo

June_Rl <- ggplot(subset(June_fluor_set_class,Class != "WT"), aes(x = Timepoint, group=Class,color = Class))+
  geom_point(aes(y = Rl), size = 2)+
  geom_line(aes(y = Rl))+ 
  geom_errorbar(aes(ymin = Rl-Rl_se, ymax = Rl+Rl_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Rl (µmol CO2 /m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))
June_Rl

June_Jo_ratio <- ggplot((subset(June_fluor_set_class, Class != "WT")), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = Jo_ratio_6800), size = 2)+
  geom_line(aes(y = Jo_ratio_6800))+ 
  geom_errorbar(aes(ymin = Jo_ratio_6800-Jo_ratio_6800_se, ymax = Jo_ratio_6800+Jo_ratio_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio (Jo/Jt)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))
June_Jo_ratio




June_photorespiration_plot <- grid.arrange(June_Jo,June_Jo_ratio, ncol=1)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_photorespiration_plot.png",plot = June_photorespiration_plot, width = 3, height = 6.5, dpi = 300)

library(cowplot)


library(RColorBrewer)
my_colors <- c("gray0","chartreuse4","indianred3","grey")
my_colors2 <- c("chartreuse4","gray0","chartreuse4","indianred3","indianred3","indianred3","gray0","gray")

library(ggsignif)
library(gridExtra)
July_assimilation <- ggplot(subset(July_summary_class, treatment == "drought" & Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = A), size = 2)+
  geom_line(aes(y = A))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se), width = 0.2)+
  ylab("Assimilation (µmol/m^2*s)")+
  xlab("Hour")+
  scale_color_manual(values = my_colors)+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_assimilation

July_gsw <- ggplot(subset(July_summary_class, treatment == "drought" & Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = gsw), size = 2)+
  geom_line(aes(y = gsw))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )
July_gsw

July_summary_class$WUE <- July_summary_class$A/July_summary_class$gsw

July_WUE <- ggplot(subset(July_summary_class, treatment == "drought"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = WUE), size = 2)+
  geom_line(aes(y = WUE, linetype = treatment))+
  scale_color_manual(values = my_colors)+
  ylab("Water use efficiency")+
  xlab("Hour")+
  theme_bw()
July_WUE


July_ETR <- ggplot(subset(July_summary_class, treatment == "drought" & Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = ETR_6800), size = 2)+
  geom_line(aes(y = ETR_6800, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR_6800-ETR_6800_se, ymax = ETR_6800+ETR_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_ETR

July_PhiPS2 <- ggplot(subset(July_summary_class,treatment == "drought" & Class != "WT"), aes(x = Timepoint, group=Class,color = Class, shape= treatment))+
  geom_point(aes(y = PhiPS2_6800,alpha=treatment), size = 2)+
  geom_line(aes(y = PhiPS2_6800,linetype = treatment,alpha=treatment))+
  scale_alpha_manual(values=c(1,0.5))+
  geom_errorbar(aes(ymin = PhiPS2_6800-PhiPS2_6800_se, ymax = PhiPS2_6800+PhiPS2_6800_se,linetype = treatment,alpha=treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("efficiency of PSII")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_PhiPS2

July_Tleaf <- ggplot(subset(July_summary_class,treatment == "drought" & Class != "WT"), aes(x = Timepoint, group=Class,color = Class, shape = treatment))+
  geom_point(aes(y = Tleaf,alpha=treatment), size = 2)+
  geom_line(aes(y = Tleaf,linetype = treatment,alpha=treatment))+ 
  scale_alpha_manual(values=c(1,0.5))+
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se,linetype = treatment,alpha=treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_Tleaf

July_PAR <- ggplot(subset(July_summary_class,treatment == "drought" & Class != "WT"), aes(x = Timepoint, group = Class, color = Class))+
  geom_point(aes(y = PAR), size = 2)+
  geom_line(aes(y = PAR))+ 
  geom_errorbar(aes(ymin = PAR-PAR_se, ymax = PAR+PAR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*sec)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

July_PAR


July_diurnal_phys_plot <- grid.arrange(July_assimilation, July_gsw,July_Tleaf,nrow=1)
July_diurnal_light_plot <- grid.arrange(July_PAR, July_PhiPS2, July_ETR ,nrow=1)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_diurnal_phys_plot.png", plot = July_diurnal_phys_plot, height = 3, width = 7, units = "in", dpi = 300)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_diurnal_light_plot.png", plot = July_diurnal_light_plot, height = 3, width = 7, units = "in", dpi = 300)


July_Jo <- ggplot(subset(July_fluor_set_class,treatment == "drought" & Class != "WT"), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_6800), size = 2)+
  geom_line(aes(y = Jo_6800, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_6800-Jo_6800_se, ymax = Jo_6800+Jo_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_Jo

July_Jo_ratio <- ggplot(subset(July_fluor_set_class,treatment == "drought" & Class != "WT"), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio_6800), size = 2)+
  geom_line(aes(y = Jo_ratio_6800, linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio_6800-Jo_ratio_6800_se, ymax = Jo_ratio_6800+Jo_ratio_6800_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

July_Jo_ratio

July_photorespiration_plot <- grid.arrange(July_Jo,July_Jo_ratio, ncol=1)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_photorespiration_plot.png",plot = July_photorespiration_plot, width = 3, height = 6.5, dpi = 300)


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

August_assimilation <- ggplot(subset(August_summary_class,Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = A,shape = treatment), size = 2)+
  geom_line(aes(y = A,linetype = treatment))+ 
  geom_errorbar(aes(ymin = A-A_se, ymax = A+A_se,linetype = treatment), width = 0.2)+
  ylab("Assimilation (µmol/m^2*s)")+
  xlab("Hour")+
  scale_color_manual(values = my_colors)+
  geom_signif(aes(x= Timepoint, y=A),comparisons = list(c("control","elite")),map_signif_level=TRUE, y_position = 15)+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

August_assimilation


August_gsw <- ggplot(subset(August_summary_class,Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = gsw,shape=treatment), size = 2)+
  geom_line(aes(y = gsw,linetype = treatment))+ 
  geom_errorbar(aes(ymin = gsw-gsw_se, ymax = gsw+gsw_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Stomatal conductance (mmol/m^2*s)")+
  xlab("Hour")+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )
August_gsw

August_ETR <- ggplot(subset(August_summary_class,Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = ETR_6800,shape=treatment), size = 2)+
  geom_line(aes(y = ETR_6800, linetype = treatment))+ 
  geom_errorbar(aes(ymin = ETR_6800-ETR_6800_se, ymax = ETR_6800+ETR_6800_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Electron transport rate (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )
August_ETR


August_PhiPS2 <- ggplot(subset(August_summary_class, Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = PhiPS2_6800,shape=treatment), size = 2)+
  geom_line(aes(y = PhiPS2_6800,linetype = treatment))+
  geom_errorbar(aes(ymin = PhiPS2_6800-PhiPS2_6800_se, ymax = PhiPS2_6800+PhiPS2_6800_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("efficiency of PSII")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

August_PhiPS2

August_Tleaf <- ggplot(subset(August_summary_class,Class != "WT"), aes(x = Timepoint,color = Class))+
  geom_point(aes(y = Tleaf), size = 2)+
  geom_line(aes(y = Tleaf,linetype = treatment))+
  geom_errorbar(aes(ymin = Tleaf-Tleaf_se, ymax = Tleaf+Tleaf_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Leaf Temp (˚C)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

August_Tleaf

August_PAR <- ggplot(subset(August_summary_class, Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = PAR), size = 2)+
  geom_line(aes(y = PAR,linetype = treatment))+ 
  geom_errorbar(aes(ymin = PAR-PAR_se, ymax = PAR+PAR_se), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("PAR (µmol/m^2*sec)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18))

August_PAR



August_diurnal_phys_plot <- grid.arrange(August_assimilation, August_gsw,August_Tleaf, nrow=1)
August_diurnal_light_plot <- grid.arrange(August_PAR, August_PhiPS2,August_ETR, nrow=1)

ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/August_diurnal_phys_plot.png", plot = August_diurnal_phys_plot, height = 3, width = 7, units = "in", dpi = 300)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/August_diurnal_light_plot.png", plot = August_diurnal_light_plot, height = 3, width = 7, units = "in", dpi = 300)

August_Jo <- ggplot(subset(August_fluor_set_class,Class != "WT"), aes(x = Timepoint, color = Class))+
  geom_point(aes(y = Jo_6800,shape=treatment), size = 2)+
  geom_line(aes(y = Jo_6800,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_6800-Jo_6800_se, ymax = Jo_6800+Jo_6800_se,linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo (µmol/m^2*s)")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )
August_Jo

August_Jo_ratio <- ggplot((subset(August_fluor_set_class,Class != "WT")), aes(x = Timepoint, color = Class, shape = treatment))+
  geom_point(aes(y = Jo_ratio_6800), size = 2)+
  geom_line(aes(y = Jo_ratio_6800,linetype = treatment))+ 
  geom_errorbar(aes(ymin = Jo_ratio_6800-Jo_ratio_6800_se, ymax = Jo_ratio_6800+Jo_ratio_6800_se, linetype = treatment), width = 0.2)+
  scale_color_manual(values = my_colors)+
  ylab("Jo ratio")+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )
August_Jo_ratio

August_photorespiration_plot <- grid.arrange(August_Jo,August_Jo_ratio, ncol=1)
ggsave(filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/August_photorespiration_plot.png",plot = August_photorespiration_plot, width = 3, height = 6.5, dpi = 300)

August_predawn_summary_class$MPa

August_predawn_summary_class$mPA
ggplot(subset(August_predawn_summary_class, Class != "WT"), aes(x= Class, y = mPA, color = treatment))+
  geom_point()+
  geom_errorbar(aes(ymin = mPA-mPA_se,ymax=mPA+mPA_se))+
  ylab("Predawn water potential (MPa)")+
  xlab(NULL)+
  theme_bw()


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

#weather dat

library(lubridate)
weather_6_27 <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
str(weather_6_27)
weather_6_27$Datetime <- ymd_hms(weather_6_27$Datetime)
weather_6_27$Date <- as.Date(weather_6_27$Datetime)
weather_6_27$Date
weather_6_27 <- subset(weather_6_27, Date == "2023-06-27")
weather_6_27$hour <- hour(weather_6_27$Datetime)


June_solar_rad <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = solar_radiation_mean), size = 2)+
  geom_line(aes(y = solar_radiation_mean))+ 
  ylab("mean solar radiation (W/m2)")+
  theme_bw()
June_solar_rad

June_air_temp <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = air_temp_mean), size = 2)+
  geom_line(aes(y = air_temp_mean))+ 
  ylab("mean temp (˚C)")+
  theme_bw()
June_air_temp

June_air_VPD <- ggplot(subset(weather_6_27, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = VPD_mean), size = 2)+
  geom_line(aes(y = VPD_mean))+ 
  ylab("VPD (kPa)")+
  theme_bw()
June_air_VPD

June_weather_station <- grid.arrange(June_solar_rad, June_air_temp, June_air_VPD, nrow=3)
ggsave(plot = June_weather_station,filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/June_weather_station.png",height = 8, width = 8, units = "in",dpi = 300 )

#July
library(lubridate)
weather_7_26 <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
str(weather_7_26)
weather_7_26$Datetime <- ymd_hms(weather_7_26$Datetime)
weather_7_26$Date <- as.Date(weather_7_26$Datetime)
weather_7_26$Date
weather_7_26 <- subset(weather_7_26, Date == "2023-07-26")
weather_7_26$hour <- hour(weather_7_26$Datetime)


July_solar_rad <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = solar_radiation_mean), size = 2)+
  geom_line(aes(y = solar_radiation_mean))+ 
  ylab("mean solar radiation (W/m2)")+
  theme_bw()
July_solar_rad

July_air_temp <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = air_temp_mean), size = 2)+
  geom_line(aes(y = air_temp_mean))+ 
  ylab("mean temp (˚C)")+
  theme_bw()
July_air_temp

July_air_VPD <- ggplot(subset(weather_7_26, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = VPD_mean), size = 2)+
  geom_line(aes(y = VPD_mean))+ 
  ylab("VPD (kPa)")+
  theme_bw()
July_air_VPD

July_weather_station <- grid.arrange(July_solar_rad, July_air_temp, July_air_VPD, nrow=3)
ggsave(plot = July_weather_station,filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/July_weather_station.png",height = 8, width = 8, units = "in",dpi = 300 )

#august

weather_8_28 <- read.csv(file = "LC_2023/2023_weather/weather_processed/2023_weather_hourly.csv")
str(weather_8_28)
weather_8_28$Datetime <- ymd_hms(weather_8_28$Datetime)
weather_8_28$Date <- as.Date(weather_8_28$Datetime)
weather_8_28$Date
weather_8_28 <- subset(weather_8_28, Date == "2023-08-28")
weather_8_28$hour <- hour(weather_8_28$Datetime)

Aug_solar_rad <- ggplot(subset(weather_8_28, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = solar_radiation_mean), size = 2)+
  geom_line(aes(y = solar_radiation_mean))+ 
  ylab("mean solar radiation (W/m2)")+
  theme_bw()
Aug_solar_rad

Aug_air_temp <- ggplot(subset(weather_8_28, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = air_temp_mean), size = 2)+
  geom_line(aes(y = air_temp_mean))+ 
  ylab("mean temp (˚C)")+
  theme_bw()
Aug_air_temp

Aug_air_VPD <- ggplot(subset(weather_8_28, hour >= 8 & hour <= 18), aes(x = hour))+
  geom_point(aes(y = VPD_mean), size = 2)+
  geom_line(aes(y = VPD_mean))+ 
  ylab("VPD (kPa)")+
  theme_bw()
Aug_air_VPD

Aug_weather_station <- grid.arrange(Aug_solar_rad, Aug_air_temp, Aug_air_VPD, nrow=3)
ggsave(plot = Aug_weather_station,filename = "LC_2023/2023_physiology_analysis/diurnal_experiments/August_weather_station.png",height = 8, width = 8, units = "in",dpi = 300 )



### summary between dates
June_fluor_set_summary_tree$treatment == "drought"
summary_total_tree <- rbind(June_fluor_set_summary_tree, July_fluor_set_summary_tree, August_fluor_set_summary_tree)

ggplot(subset(summary_total_tree, Ci > 1 & Ci < 500 & Jo_ratio > 0 & Jo_ratio < 1.0), aes(x=Ci,y=Jo_ratio,color=Timepoint))+
  geom_point()

ggplot(subset(summary_total_tree, Jo_ratio > 0 & Jo_ratio < 1.0), aes(x=gsw,y=Jo_ratio,color=Timepoint))+
  geom_point()



June_fluor_set_class$sample_date <- "06-27-2023"
June_fluor_set_class$treatment <- "drought"
June_fluor_set_class_II <- subset(June_fluor_set_class, select = c(Class,treatment,n,Timepoint,A,A_se,PAR,PAR_se,Rd,Rd_se,Tleaf,Tleaf_se,gsw,gsw_se,ETR_6800,ETR_6800_se,PhiPS2_6800,PhiPS2_6800_se,Jo_6800,Jo_6800_se,Jo_ratio_6800,Jo_ratio_6800_se,Ci,Ci_se,sample_date))

July_fluor_set_class$sample_date <- "07-26-2023"
July_fluor_set_class_II <- subset(July_fluor_set_class, select = c(Class,treatment,n,Timepoint,A,A_se,PAR,PAR_se,Rd,Rd_se,Tleaf,Tleaf_se,gsw,gsw_se,ETR_6800,ETR_6800_se,PhiPS2_6800,PhiPS2_6800_se,Jo_6800,Jo_6800_se,Jo_ratio_6800,Jo_ratio_6800_se,Ci,Ci_se,sample_date))

August_fluor_set_class$sample_date <- "08-28-2023" 
August_fluor_set_class_II <- subset(August_fluor_set_class, select = c(Class,treatment,n,Timepoint,A,A_se,PAR,PAR_se,Rd,Rd_se,Tleaf,Tleaf_se,gsw,gsw_se,ETR_6800,ETR_6800_se,PhiPS2_6800,PhiPS2_6800_se,Jo_6800,Jo_6800_se,Jo_ratio_6800,Jo_ratio_6800_se,Ci,Ci_se,sample_date))


total_diurnal_dat <- rbind(June_fluor_set_class_II,July_fluor_set_class_II,August_fluor_set_class_II)
#rename some columns


total_diurnal_dat <- rename(total_diurnal_dat, ETR6800 = ETR_6800, ETR6800_se = ETR_6800_se, PhiPS26800 = PhiPS2_6800, PhiPS26800_se = PhiPS2_6800_se,Jo6800 = Jo_6800,Jo6800_se=Jo_6800_se,Joratio6800 = Jo_ratio_6800,Joratio6800_se = Jo_ratio_6800_se)

#convert to longer
total_diurnal_dat_longer_1 <- pivot_longer(total_diurnal_dat, names_to = "meas",values_to = "value", cols = c(A,PAR,Rd,Tleaf,gsw,ETR6800,PhiPS26800,Jo6800,Joratio6800))
total_diurnal_dat_longer_1 <-subset(total_diurnal_dat_longer_1,select = c("Class","treatment","n","sample_date","Timepoint","meas","value"))

total_diurnal_dat_longer_2 <- pivot_longer(total_diurnal_dat, names_to = "meas",values_to = "se", cols = c(A_se,PAR_se,Rd_se,Tleaf_se,gsw_se,ETR6800_se,PhiPS26800_se,Jo6800_se,Joratio6800_se))
total_diurnal_dat_longer_2 <-subset(total_diurnal_dat_longer_2,select = c("Class","treatment","n","sample_date","Timepoint","meas","se"))
library(stringr)
total_diurnal_dat_longer_2$meas <- str_replace(total_diurnal_dat_longer_2$meas, "_se", "")


diurnal_longer_clean <- inner_join(total_diurnal_dat_longer_1,total_diurnal_dat_longer_2, by = c("Class","treatment","sample_date","n","Timepoint","meas"))

diurnal_longer_clean2 <- diurnal_longer_clean %>% distinct(value, .keep_all = TRUE)

#t <- pivot_longer(total_diurnal_dat, names_to = c("measurement"),names_sep = "_",values_to= "value",cols = c(A,A_se,PAR,PAR_se,Rd,Rd_se,Tleaf,Tleaf_se,gsw,gsw_se,ETR_6800,ETR_6800_se,PhiPS2_6800,PhiPS2_6800_se,Jo_6800,Jo_6800_se,Jo_ratio_6800,Jo_ratio_6800_se,Ci,Ci_se))


ggplot(subset(diurnal_longer_clean,treatment =="drought"), aes(x=Timepoint,y=value,color=sample_date,shape = Class))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = my_colors)+
  facet_wrap(~meas, scales = "free")

ggplot(subset(diurnal_longer_clean,treatment =="drought"), aes(x=Timepoint,y=value,color=Class,shape = sample_date))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = my_colors)+
  facet_wrap(sample_date~meas, scales = "free")






#summarizing over sample dates
total_diurnal_summary <- total_diurnal_dat %>% group_by(Timepoint,Class,treatment) %>% summarize(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = gsw_sd/(sqrt(n)),
  ETR_sd = sd(ETR6800),
  ETR = mean(ETR6800),
  ETR_se = ETR_sd/(sqrt(n)),
  PhiPS2_sd = sd(PhiPS26800),
  PhiPS2 = mean(PhiPS26800),
  PhiPS2_se = PhiPS2_sd/(sqrt(n)),
  PAR_sd = sd(PAR),
  PAR = mean(PAR),
  PAR_se = PAR_sd/(sqrt(n)),
  Jo_sd = sd(Jo6800),
  Jo = mean(Jo6800),
  Jo_se = Jo_sd/(sqrt(n)),
  Joratio_sd = sd(Joratio6800),
  Joratio = mean(Joratio6800),
  Joratio_se = Joratio_sd/(sqrt(n)),
  Tleaf_sd = sd(Tleaf),
  Tleaf = mean(Tleaf),
  Tleaf_se = Tleaf_sd/(sqrt(n))
)

#converting to long

diurnal_summ_longer_1 <- pivot_longer(total_diurnal_summary, names_to = "meas",values_to = "value", cols = c(A,PAR,Tleaf,gsw,ETR,PhiPS2,Jo,Joratio))
diurnal_summ_longer_1 <-subset(diurnal_summ_longer_1,select = c("Class","treatment","n","Timepoint","meas","value"))

diurnal_summ_longer_2 <- pivot_longer(total_diurnal_summary, names_to = "meas",values_to = "se", cols = c(A_se,PAR_se,Tleaf_se,gsw_se,ETR_se,PhiPS2_se,Jo_se,Joratio_se))
diurnal_summ_longer_2 <-subset(diurnal_summ_longer_2,select = c("Class","treatment","n","Timepoint","meas","se"))

diurnal_summ_longer_2$meas <- str_replace(diurnal_summ_longer_2$meas, "_se", "")

diurnal_summ_longer <- inner_join(diurnal_summ_longer_1,diurnal_summ_longer_2, by = c("Class","treatment","n","Timepoint","meas"))

diurnal_summ_longer$meas <- factor(diurnal_summ_longer$meas, levels = c("A","gsw","Tleaf","PAR","PhiPS2","ETR","Jo","Joratio"))

diurnal_summary_plot <- ggplot(subset(diurnal_summ_longer,treatment =="drought" & Class != "WT"), aes(x=Timepoint,y=value,color=Class))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=value-se,ymax=value+se),width=0.2)+
  scale_color_manual(values = my_colors)+
  facet_wrap(~meas, scales = "free",nrow = 4)+
  xlab("Hour")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(
    breaks = seq(8,18,by=2),
    limits = c(8,18)
  )

ggsave(filename="LC_2023/2023_physiology_analysis/diurnal_experiments/diurnal_summary_plot.png",plot=diurnal_summary_plot, height = 8, width = 12, units = "in", dpi = 300)
  


  