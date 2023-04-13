##Metabolics 2023
library(ggplot2)
library(gridExtra)
library((ggthemes))
library(dplyr)
library(ggsignif)
library(PCAtools)
library(GGally)
library(ggbiplot)


str(metabolomics)
metabolomics <- read.csv(file = "2022_molecular_analysis/2022_metabolomics/metabolomics_final_raw.csv")
str(metabolomics)
growth <- read.csv(file = "2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
metabolomics_II <- inner_join(metabolomics, growth, by = "ID")



metabolomics_III <- inner_join(metabolomics_II, by_tree, by = "ID")
str(metabolomics_III)  
##looking at data
glycolate_event <- ggplot(metabolomics_III, aes(x=event_short, y = glycolic.acid, fill = construct2.x))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("13-15E","CT3")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("2H","CT3")),map_signif_level=TRUE, y_position = 6000)

glycolate_event
ggsave(filename = "glycolate_event.png", plot = glycolate_event, dpi = 300)


#gene expression data from both years
#read in compiled csv
by_tree <- read.csv("PAG_analysis/Gene_exp_bytree_both_yrs.csv")
str(by_tree)
by_tree$PLGG1_2021 <- as.numeric(by_tree$PLGG1_2021)
by_tree$GDH_2021 <- as.numeric(by_tree$GDH_2021)
by_tree$MS_2021 <- as.numeric(by_tree$MS_2021)
by_tree$PLGG1_avg <- as.numeric(by_tree$PLGG1_avg)
by_tree$GDH_avg <- as.numeric(by_tree$GDH_avg)
by_tree$MS_avg <- as.numeric(by_tree$MS_avg)

metabolomics_IV <- inner_join(metabolomics_III, by_tree, by = "ID")

glycolate_plgg1_plot <- ggplot(metabolomics_IV,aes(x= PLGG1_avg.y, y = log(glycolic.acid)))+
  geom_point(aes(colour = metabolomics_IV$event_short, shape = construct2),size=2)
 
glycolate_plgg1_lm <- lm(log(glycolic.acid)~PLGG1_avg.y + event_short, data = metabolomics_IV)
glycolate_plgg1_aov <- aov(log(glycolic.acid)~PLGG1_avg.y + event_short, data = metabolomics_IV)
plot(glycolate_plgg1_aov, which=1)

anova(glycolate_plgg1_lm)


glycolate_plgg1_plot
ggsave(filename = "glycolate_plgg1_plot.png", plot = glycolate_plgg1_plot, dpi = 300)


ggplot(metabolomics_IV, aes(x=V419, y=glycolic.acid))+
  geom_point(aes(colour = event_short))

####PCA of metabolites ########

#PCA_set <- metabolomics_II[,2:74]


PCA_set <- na.omit(metabolomics)
PCA_set <- PCA_set[,3:67]

PCA_metabolomics <- prcomp(PCA_set, center = TRUE, scale. = TRUE)

ggbiplot(PCA_metabolomics)


PCA_set_II <- metabolomics[,c(14,43,50,36,55,67,48,38,15,53,21,35,8,4)]
PCA_labels <- metabolomics_II[,c(1,14,43,50,36,55,67,48,38,15,53,21,35,8,4,71,72)]
PCA_set_II <- na.omit(PCA_set_II)
PCA_labels <- na.omit(PCA_labels)

PCA_metabolomics_II <- prcomp(PCA_set_II, center = TRUE, scale. = TRUE)
ggbiplot(PCA_metabolomics_II, ellipse =TRUE, labels = rownames(PCA_labels), groups = PCA_labels$event_short)


ggbiplot(PCA_tree2, choices = c(2,3),obs.scale = 1,var.scale = 1,ellipse = TRUE, labels=rownames(transgene_by_tree), groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")

phys <- read.csv("../LC_June2022/Li6800_data/Aci_parameters_list.csv")

growth <- read.csv("../LC_Nov_22_update/DBH_H_timeline_CT1_excluded_9_22.csv")

str(metabolomics)
metabolomics_II <- inner_join(metabolomics, growth, by = "ID")

metabolomics_III <- inner_join(metabolomics_II, by_tree, by = "ID")

##looking at data, expected relationships
ggplot(metabolomics_III, aes(x=event_short, y = glycolic.acid))+
  geom_boxplot()

ggplot(metabolomics_III, aes(x=event_short, y = glyceric.acid))+
  geom_boxplot()

ggplot(metabolomics_III, aes(x=event_short, y = L.serine))+
  geom_boxplot()

ggplot(metabolomics_III, aes(x=event_short, y = D.malic.acid))+
  geom_boxplot()

ggplot(metabolomics_III, aes(x=event_short, y = glycine))+
  geom_boxplot()


##########



phys$ID <- phys$tree
physII <- read.csv("../LC_June2022/Li6800_data/Aci_parameters_list_flr.csv")
physII$ID <- physII$tree
growth <- read.csv("../LC_Nov_22_update/DBH_H_timeline_CT1_excluded_9_22.csv")

phys <- left_join(phys, growth, by = "ID")
physII <- left_join(physII, growth, by = "ID")

#find event means for phys and then growth variables
phys_event <- phys %>% group_by(event)
phys_summary <- phys_event %>% dplyr::summarise(
  n = n(),
  Ci_sd = sd(Ci_A_410),
  Ci = mean(Ci_A_410),
  Ci_se = (Ci_sd/n),
  Amax_sd = sd(A_A_410),
  Amax = mean(A_A_410),
  Amax_se = (Amax_sd/n),
  Vcmax_sd = sd(Vcmax),
  Vcmax = mean(Vcmax),
  Vcmax_se = (Vcmax_sd/n),
  Jmax_sd = sd(Jmax),
  Jmax = mean(Jmax),
  Jmax_se = (Jmax_sd/n),
  Jmax_Vcmax_ratio_sd = sd(Jmax_Vcmax_ratio),
  Jmax_Vcmax = mean(Jmax_Vcmax_ratio),
  Jmax_Vcmax_ratio_se = (Jmax_Vcmax_ratio_sd/n),
  Rd_sd = sd(Rd),
  Rd = mean(Rd),
  Rd_se = (Rd_sd/n)
)

phys_event_II <- physII %>% group_by(event)


physII_summary <- phys_event_II %>% dplyr::summarise(
  n = n(),
  PhiPS2_PhiCO2_sd = sd(PhiPS2_PhiCO2_A_410),
  PhiPS2_PhiCO2 = mean(PhiPS2_PhiCO2_A_410),
  PhiPS2_PhiCO2_se = (PhiPS2_PhiCO2_sd/n)
)

phys_summary_III <- inner_join(phys_summary, physII_summary, by = "event")

growth_event <- growth %>% group_by(event)

growth_summary <- growth_event %>% dplyr::summarize(
  n = n(),
  VI_sd = sd(V419),
  VI = mean(V419),
  VI_se = (VI_sd/n),
  Chl_sd = sd(Aug_SPAD),
  Chl = mean(Aug_SPAD),
  Chl_se = (Chl_sd/n)
)

#join datasets of event means

growth_phys_summary <- inner_join(phys_summary_III, growth_summary, by = "event")

metabolomics_summary <- inner_join(metabolomics, growth_phys_summary, )