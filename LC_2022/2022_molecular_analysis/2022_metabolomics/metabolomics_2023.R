##Metabolics 2023
library(ggplot2)
library(gridExtra)
library((ggthemes))
library(dplyr)
library(ggsignif)
library(PCAtools)
library(GGally)
library(ggbiplot)



metabolomics <- read.csv(file = "2022_molecular_analysis/2022_metabolomics/metabolomics_final_raw.csv")
str(metabolomics)
metabolomics$ribitol..IS. <- as.numeric(metabolomics$ribitol..IS.)
growth <- read.csv(file = "2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
metabolomics_II <- inner_join(metabolomics, growth, by = "ID")


#gene expression data by tree
#gene expression data from both years
#read in compiled csv
by_tree <- read.csv("2022_molecular_analysis/2022_gene_expression/Gene_exp_bytree_both_yrs.csv")
str(by_tree)
by_tree$PLGG1_2021 <- as.numeric(by_tree$PLGG1_2021)
by_tree$GDH_2021 <- as.numeric(by_tree$GDH_2021)
by_tree$MS_2021 <- as.numeric(by_tree$MS_2021)
by_tree$PLGG1_avg <- as.numeric(by_tree$PLGG1_avg)
by_tree$GDH_avg <- as.numeric(by_tree$GDH_avg)
by_tree$MS_avg <- as.numeric(by_tree$MS_avg)
by_tree$transgene_avg <- as.numeric(by_tree$transgene_avg)



metabolomics_III <- inner_join(metabolomics_II, by_tree, by = "ID")
str(metabolomics_III)  
##looking at data
glycolate_event <- ggplot(metabolomics_III, aes(x=event_short, y = glycolic.acid, fill = construct2.x))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("13-15E","CT3")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("2H","CT3")),map_signif_level=TRUE, y_position = 6000)+
  ylab("Relative glycolate level (to IS)")+
  xlab("event")

glycolate_event
ggsave(filename = "glycolate_event.png", plot = glycolate_event, dpi = 300)

###statistics for glycolate#####

##event comparison###

# Combine 16-20 and 8-9D
metabolomics_III$event2 <- metabolomics_III$event
metabolomics_III$event2[metabolomics_III$event == "LC-102 16-20" |
                metabolomics_III$event == "LC-102 8-9D"] <- "escape"

glycolate_event_model <- lm(glycolic.acid ~ event2 +block, data = metabolomics_III)
summary(glycolate_event_model)
#block is significant, include as random effect

library(nlme)
library(emmeans)

glycolate_event_model.2 <- lme(log(glycolic.acid) ~ event2, random = ~1|block, data = metabolomics_III)
summary(glycolate_event_model.2)

plot(fitted(glycolate_event_model.2), residuals(glycolate_event_model.2), xlab="Fitted Values",
     ylab="Studentized Residuals",
     main="Fitted vs. Residuals"); abline(h=0)
qqnorm(residuals(glycolate_event_model.2)); qqline(residuals(glycolate_event_model.2))
##log transformation definitely helped

emmeans(glycolate_event_model.2, specs = pairwise ~event2)

emmeans(glycolate_event_model.2, specs = pairwise ~event2)$contrasts


glycolate_event_summary <- as.data.frame(emmeans(glycolate_event_model.2, specs = pairwise ~event2)$contrasts)
glycolate_event_summary$effect_size <- (round(exp(-1*(glycolate_event_summary$estimate)),3)-1)*100
glycolate_event_summary$p.value <- round(glycolate_event_summary$p.value,3)
glycolate_event_summary$estimate <- round(glycolate_event_summary$estimate,3)
glycolate_event_summary$SE <- round(glycolate_event_summary$SE,2)



write.csv(glycolate_event_summary, file = "2022_molecular_analysis/2022_metabolomics/glycolate_event_summary_stats.csv")

metabolomics_IV <- inner_join(metabolomics_III, by_tree, by = "ID")

glycolate_plgg1_plot <- ggplot(metabolomics_IV,aes(x= PLGG1_avg.y, y = log(glycolic.acid)))+
  geom_point(aes(colour = metabolomics_III$event_short, shape = construct2),size=2)
 
glycolate_plgg1_lm <- lm(log(glycolic.acid)~PLGG1_avg.y + event_short, data = metabolomics_IV)
glycolate_plgg1_aov <- aov(log(glycolic.acid)~PLGG1_avg.y + event_short, data = metabolomics_IV)
plot(glycolate_plgg1_aov, which=1)

anova(glycolate_plgg1_lm)


glycolate_plgg1_plot
ggsave(filename = "glycolate_plgg1_plot.png", plot = glycolate_plgg1_plot, dpi = 300)


ggplot(metabolomics_IV, aes(x=V419, y=glycolic.acid))+
  geom_point(aes(colour = event_short))

####PCA of metabolites ########

##part of purpose is to explore feature reduction

PCA_metabolites_set <- metabolomics_III[,c(3:67,73)]
write.csv(PCA_metabolites_set,file = "Metabolites_raw_event.csv")
#only looking at metabolites that were real numbers for all data
#should ask Connor or Jeff if NAs can be interpreted as 0s
PCA_metabolites_set <- PCA_metabolites_set %>%
  select_if(~ !any(is.na(.)))

#exporting CsV for python
PCA_metabolites_set_II <- metabolomics_III[,c(2:67,73)]
PCA_metabolites_set_II <- PCA_metabolites_set_II %>%
  select_if(~ !any(is.na(.)))

PCA_meta <- subset(metabolomics_III, select = c("ID","event_short","block","PLGG1_2022","GDH_2022","MS_2022"))
PCA_metabolites_export <- inner_join(PCA_metabolites_set_II, PCA_meta, by = "ID")
write.csv(PCA_metabolites_export, file = "PCA_metabolites_export.csv")
##


PCA_metabolites <- prcomp(PCA_metabolites_set[,c(1:53)], center = TRUE, scale. = TRUE, print = TRUE)


#which PCs are most important?
ggscreeplot(PCA_metabolites)

ggbiplot(PCA_metabolites)
#The main distinguisher is that some trees seem to have less of everything
#sucrose is set off to the side as well
ggbiplot(PCA_metabolites, choices = c(1,3))
ggbiplot(PCA_metabolites, choices = c(1,4))
ggbiplot(PCA_metabolites, choices = c(2,3))
ggbiplot(PCA_metabolites, choices = c(2,4))

##how do different events compare in PC1 versus PC2 space
ggbiplot(PCA_metabolites,obs.scale = 1, var.axes = FALSE, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_metabolites_set$event_short )+
  ggtitle("Principal components of poplar leaf metabolome")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (42.1%)") + ylab ("PC2 (9.4%)")+
  geom_point(aes(colour= PCA_metabolites_set$event_short), size = 1)

ggsave(filename = "2022_molecular_analysis/2022_metabolomics/metabolite_PCA.png")
#no one event is clustered away from the rest. 

PCA_2_3 <- ggbiplot(PCA_metabolites,obs.scale = 1, choices = c(2,3), var.axes = FALSE, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_metabolites_set$event_short )+
  ggtitle("Principal components of poplar leaf metabolome")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC2 (9.4%)") + ylab("PC3 (6.5%)")
  geom_point(aes(colour= PCA_metabolites_set$event_short), size = 1)

ggsave(filename = "2022_molecular_analysis/2022_metabolomics/metabolite_PCA_2_3.png", dpi = 300)

## Use PCA loadings to infer metabolites contributing to PCs

PCA_loadings <- as.data.frame(PCA_metabolites$rotation)
PCA_loadings$metabolite <- row.names(PCA_loadings)
#remove internal standard ribose
PCA_loadings <- subset(PCA_loadings, metabolite != "ribitol..IS.")

#PC1
PC1_max <- slice_max(PCA_loadings,order_by = PC1, n= 10)
PC1_max <- subset(PC1_max, select = c("metabolite","PC1"))
PC1_min <- slice_min(PCA_loadings,order_by = PC1, n=10)
PC1_min <- subset(PC1_min, select = c("metabolite","PC1"))

PC1_loadings <- rbind(PC1_max, PC1_min)

PC1_min_graph <- ggplot(subset(PC1_loadings, PC1 < -0.1), aes(x=reorder(metabolite, PC1), y = PC1))+
  geom_bar(stat='identity', fill = "orange3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC1 loading")+
  xlab("metabolite")

PC1_min_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC1_loadings_min.png",plot = PC1_min_graph,  height = 5, width = 4, units = "in", dpi = 300)

PC1_max_graph <- ggplot(subset(PC1_loadings, PC1 > -0.1), aes(x=reorder(metabolite, -1*PC1), y = PC1))+
  geom_bar(stat='identity', fill = "blue3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC1 loading")+
  xlab("metabolite")

PC1_max_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC1_loadings_max.png",plot = PC1_max_graph, height = 5, width = 4, units = "in", dpi = 300)

##PC2

PC2_max <- slice_max(PCA_loadings,order_by = PC2, n= 10)
PC2_max <- subset(PC2_max, select = c("metabolite","PC2"))
PC2_min <- slice_min(PCA_loadings,order_by = PC2, n=10)
PC2_min <- subset(PC2_min, select = c("metabolite","PC2"))

PC2_loadings <- rbind(PC2_max, PC2_min)

PC2_min_graph <- ggplot(subset(PC2_loadings, PC2 < 0), aes(x=reorder(metabolite, PC2), y = PC2))+
  geom_bar(stat='identity', fill = "orange3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC2 loading")+
  xlab("metabolite")

PC2_min_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC2_loadings_min.png",plot = PC2_min_graph,  height = 5, width = 4, units = "in", dpi = 300)

PC2_max_graph <- ggplot(subset(PC2_loadings, PC2 > 0), aes(x=reorder(metabolite, -1*PC2), y = PC2))+
  geom_bar(stat='identity', fill = "blue3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC2 loading")+
  xlab("metabolite")

PC2_max_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC2_loadings_max.png",plot = PC2_max_graph, height = 5, width = 4, units = "in", dpi = 300)


#PC3
PC3_max <- slice_max(PCA_loadings,order_by = PC3, n= 10)
PC3_max <- subset(PC3_max, select = c("metabolite","PC3"))
PC3_min <- slice_min(PCA_loadings,order_by = PC3, n=10)
PC3_min <- subset(PC3_min, select = c("metabolite","PC3"))

PC3_loadings <- rbind(PC3_max, PC3_min)

PC3_min_graph <- ggplot(subset(PC3_loadings, PC3 < 0), aes(x=reorder(metabolite, PC3), y = PC3))+
  geom_bar(stat='identity', fill = "orange3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC3 loading")+
  xlab("metabolite")

PC3_min_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC3_loadings_min.png",plot = PC3_min_graph,  height = 5, width = 4, units = "in", dpi = 300)


PC3_max_graph <- ggplot(subset(PC3_loadings, PC3 > 0), aes(x=reorder(metabolite, -1*PC3), y = PC3))+
  geom_bar(stat='identity', fill = "blue3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC3 loading")+
  xlab("metabolite")


PC3_max_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC3_loadings_max.png",plot = PC3_max_graph, height = 5, width = 4, units = "in", dpi = 300)




##assign event level mean and se of PC1 and PC2
metabolites_PCA_score <- as.data.frame(PCA_metabolites$x)

metabolites_PCA_score$ID <- PCA_metabolites_set$ID
metabolites_PCA_score$event_short <- PCA_metabolites_set$event_short

metabolites_PCA_score_event <- metabolites_PCA_score %>% group_by(event_short)

metabolites_event_summary <- metabolites_PCA_score_event %>% dplyr::summarise(
  n = n(),
  PC_Exp_1_sd = sd(PC1),
  PC_Exp_1 = mean(PC1),
  PC_Exp_1_se = (PC_Exp_1_sd/(sqrt(n))),
  PC_Exp_2_sd = sd(PC2),
  PC_Exp_2 = mean(PC2),
  PC_Exp_2_se = (PC_Exp_2_sd/(sqrt(n))))


PCA_rotations <- as.data.frame(PCA_metabolites$rotation)



#at individual tree level looking at metabolites and gene expression

PCA_set <- na.omit(metabolomics_III)

PCA_set_all <- PCA_set[,c(3:67,100,120,121,122)]

PCA_metabolomics <- prcomp(PCA_set_all, center = TRUE, scale. = TRUE)

ggbiplot(PCA_metabolomics)

##Plot with ellipses around events

ggbiplot(PCA_metabolomics, obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_set$event_short)+
  ggtitle("PCA of molecular traits by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (34.4%)") + ylab ("PC2 (11.4%)")+
  geom_point(aes(colour= PCA_set$event_short), size = 1)


PCA_set_transgenic <- subset(PCA_set, event_short == "13-15B" | event_short == "13-15E" | event_short == "1C" | event_short == "5A" | event_short == "2H" | event_short == "5C" | event_short == "4A" | event_short == "7")

PCA_set_transgenic_id <- PCA_set_transgenic[,c(68:73)]

PCA_set_transgenic <- PCA_set_transgenic[,c(3:67,100,120,121,122)]

PCA_metabolomics_transgenic <- prcomp(PCA_set_transgenic, center = TRUE, scale. = TRUE)

ggbiplot(PCA_metabolomics_transgenic)
ggbiplot(PCA_metabolomics_transgenic, obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_set_transgenic_id$event_short)+
  ggtitle("PCA of molecular traits by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (31.7%)") + ylab ("PC2 (13.9%)")+
  geom_point(aes(colour= PCA_set_transgenic_id$event_short), size = 1)


####event level comparison between transgenics

PCA_set_event <- PCA_set[,c(3:67,73,100,121,122)]

PCA_set_event <- PCA_set_event %>% group_by(event_short)
#remove controls and events with sample size less than 2


PCA_event_summary <- dplyr::summarise_all(PCA_set_event, list(mean =mean,n = length,sd = sd))
#remove controls and events with sample size less than 2
PCA_event_summary_transgenic <- PCA_event_summary[c(-1,-3,-10,-11),] 

PCA_metabolomics_event <- prcomp(PCA_event_summary_transgenic[,-1], center = TRUE, scale. = TRUE)

ggbiplot(PCA_metabolomics_event)

ggbiplot(PCA_metabolomics_event, obs.scale = 20, var.scale = 14, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_event_summary_transgenic$event_short)+
  ggtitle("PCA of molecular traits by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (50.5%)") + ylab ("PC2 (16.2%)")+
  geom_point(aes(colour= PCA_event_summary_transgenic$event_short), size = 1)





#looking at metabolites that were significantly different between transgenic and control
PCA_set_II <- metabolomics[,c(14,43,50,36,55,67,48,38,15,53,21,35,8,4)]
PCA_labels <- metabolomics_II[,c(1,14,43,50,36,55,67,48,38,15,53,21,35,8,4,71,72)]
PCA_set_II <- na.omit(PCA_set_II)
PCA_labels <- na.omit(PCA_labels)

PCA_metabolomics_II <- prcomp(PCA_set_II, center = TRUE, scale. = TRUE)
ggbiplot(PCA_metabolomics_II, ellipse =TRUE, groups = PCA_labels$event_short)



feature_importances <- read.csv("2022_molecular_analysis/2022_metabolomics/metab_feat_importance.csv")


str(feature_importances)

metabolite_importances <- ggplot(feature_importances, aes(y = MDI, x = reorder(Name,-MDI)))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Metabolite")+
  ylab("Impact on model (mean decrease in impurity)")+
  theme(text = element_text(size=20))
  
metabolite_importances

ggsave(filename = "metabolites_importances.png", plot = metabolite_importances, width = 16, height = 8, units = "in", dpi = 300)


  