## analysis for PAG

#first step - perform PCA w/ tree level data to generate PCs that capture main components of variation


library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
library(PCAtools)
library(GGally)


#Start by looking at the gene expression data from both years
#read in compiled csv
by_tree <- read.csv("Gene_exp_bytree_both_yrs.csv")
str(by_tree)
by_tree$PLGG1_2021 <- as.numeric(by_tree$PLGG1_2021)
by_tree$GDH_2021 <- as.numeric(by_tree$GDH_2021)
by_tree$MS_2021 <- as.numeric(by_tree$MS_2021)
by_tree$PLGG1_avg <- as.numeric(by_tree$PLGG1_avg)
by_tree$GDH_avg <- as.numeric(by_tree$GDH_avg)
by_tree$MS_avg <- as.numeric(by_tree$MS_avg)

#unfortunately only have data at the tree level in both years for a subset of the large block
#Trees in small block will not have any PLGG1 data from 2021 - omit them
by_tree_clean <- subset(by_tree, PLGG1_2021 > 0)
#For much of the analysis, only going to look at relationships amongst transgenic trees
transgene_by_tree <- subset(by_tree_clean, construct2 == "transgenic")

#compare expression for each gene from year to year - to see how consistent generally
ggplot(by_tree_clean, aes(x=PLGG1_2021, y=PLGG1_2022, color=Event))+
  geom_point()
ggplot(by_tree_clean, aes(x= GDH_2021, y=GDH_2022, color = Event))+
  geom_point()
ggplot(by_tree_clean, aes(x= MS_2021, y=MS_2022, color = Event))+
  geom_point()
#MS seems more scattered than GDH

#What is the level of correlation
tree_level_corr_matrix <- (by_tree_clean[,6:11])
str(tree_level_corr_matrix)
ggpairs(tree_level_corr_matrix)
#PLGG1 - 0.699
#GDH - 0.772
#MS - 0.761

#Overall best to include all data in PCA.



#run PCA on gene expression data inputting both years


PCA_treeset <- transgene_by_tree[,6:11]
PCA_tree <- prcomp(PCA_treeset, center = TRUE, scale. = TRUE)
summary(PCA_tree)
#plot w/ vectors
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

#detach("package:reshape2", unload=TRUE)
#detach("package:reshape", unload=TRUE)

ggbiplot(PCA_tree)
#cool, PLGG1 seems to be running in different direction than the two transgenes. But where are the events here?
#with IDs and event clusters
tree_PCA <- ggbiplot(PCA_tree, obs.scale = 20, var.scale = 14, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (35.0%)") + ylab ("PC2 (32.1%)")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)
  

tree_PCA
ggsave(filename = "tree_PCA.png", plot = tree_PCA, width = 8, height = 4, units = "in",dpi = 300)

#2H, 13-15E cluster together, PC1 seems to be defined by PLLG1 expression. 5C is squarely in the middle.

tree_PCA2 <- ggbiplot(PCA_tree, choices = c(3,4), obs.scale = 1.1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC3 (18.0%)")+
  ylab("PC4 (7.8%)")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

tree_PCA2
ggsave(filename = "tree_PCA2.png", plot = tree_PCA2, width = 8, height = 4, units = "in", dpi = 300)

tree_PCA3 <- ggbiplot(PCA_tree, choices = c(2,3), obs.scale = 1.1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

tree_PCA3
ggsave(filename = "tree_PCA3.png", plot = tree_PCA3, width = 8, height = 4, units = "in", dpi = 300)


tree_PCA4 <- ggbiplot(PCA_tree, choices = c(1,4), obs.scale = 1.1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

tree_PCA4
ggsave(filename = "tree_PCA4.png", plot = tree_PCA4, width = 8, height = 4, units = "in",dpi = 300)
#Cool, PCs seem useful in describing the variation of the data more efficiently. They also suggest some differences between events
#But which PCs to include in downstream analysis to summarize gene expression data?
#https://rpubs.com/Bury/ClusteringOnPcaResults suggests a method. Choosing PCs with eigenvalues greater to or equal than 1
#Another site - some forum digging suggested that the sdev output squared very closely approximates the eigenvalues
#Think this is ok as R documentation states:
#
#sdev = the standard deviations of the principal components (i.e., the square roots of the eigenvalues of the covariance/correlation matrix, though the calculation is actually done with the singular values of the data matrix).

PCA_tree$sdev^2
# the first three are greater than 1. This is what I intuitively would have guessed as well

#Now need to assign the rotation values to individual trees and then summarize event means
PCA_score <- as.data.frame(PCA_tree$x)

PCA_score$ID <- transgene_by_tree$ID
PCA_score$Event <- transgene_by_tree$Event

PCA_score_event <- PCA_score %>% group_by(Event)

gene_summary <- PCA_score_event %>% dplyr::summarise(
  n = n(),
  PC1_sd = sd(PC1),
  PC1 = mean(PC1),
  PC1_se = (PC1_sd/n),
  PC2_sd = sd(PC2),
  PC2 = mean(PC2),
  PC2_se = (PC2_sd/n),
  PC3_sd = sd(PC3),
  PC3 = mean(PC3),
  PC3_se = (PC3_sd/n)
      
)

#bring in physiology and growth data

phys <- read.csv("../LC_June2022/Li6800_data/Aci_parameters_list.csv")
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

gene_summary$event <- gene_summary$Event
full_summary <- inner_join(gene_summary, growth_phys_summary, by = "event")

#
#examine relationship between variables via correlation matrix

full_corrset <- full_summary[,c(4,7,10,15,18,21,24,27,30,34,38,41)]
full_cor_plot <- ggpairs(full_corrset, upper = list(continuous = wrap("cor", size = 2)), lower = list(continuous = wrap("points", size = 0.5)), diag = list(continuous ="blankDiag"), axisLabels = "none")
ggsave(filename = "full_cor_plot.png", plot = full_cor_plot, width = 6, height = 4, dpi = 300)
#ok, pretty busy


full_set <- full_summary[,c(4,7,10,15,18,21,24,27,30,34,38,41)]
PCA_full <- prcomp(full_set, center = TRUE, scale. = TRUE)
summary(PCA_full)

ggbiplot(PCA_full)
#nice

full_summary$Event_short <- c("13-15B","13-15E","2H","4A","5A","5C","7")


library(RColorBrewer)
myColors <- c('red',
              'orange',
              'green2',
              'green4',
              'blue',
              'purple',
              'pink'
)


names(myColors) <- levels(full_summary$Event_short)
colorScale <- scale_fill_manual(name = "Event_short",values = myColors)


full_PCA <- ggbiplot(PCA_full, obs.scale = 1.5, var.scale = 6.5, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC1 (40.5%)")+
  ylab("PC2 (28.4%")+
  geom_point(aes(colour=full_summary$Event_short), size = 3)

full_PCA 
ggsave(filename = "full_PCA.png", plot = full_PCA, width = 8, height = 4, units = "in", dpi = 300)


full_PCA2 <- ggbiplot(PCA_full, choices = c(3,4), obs.scale = 1, var.scale = 1.6, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=full_summary$Event_short), size = 3)+
  xlab("PC3 (16.6%)")+
  ylab("PC4 (8.0%)")
  
 

full_PCA2
ggsave(filename = "full_PCA2.png", plot = full_PCA2, width = 8, height = 4, units = "in", dpi = 300)

full_PCA3 <- ggbiplot(PCA_full, choices = c(2,3), obs.scale = 1, var.scale = 3, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=full_summary$Event_short), size = 3)+
  xlab("PC2 (28.3%)")+
  ylab("PC3 (16.6%)")

full_PCA3
ggsave(filename = "full_PCA3.png", plot = full_PCA3, width = 8, height = 4, units = "in", dpi = 300)

full_PCA4 <- ggbiplot(PCA_full, choices = c(1,4), obs.scale = 1, var.scale = 1.5, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=full_summary$Event_short), size = 3)+
  xlab("PC1 (40.5%)")+
  ylab("PC4 (8.0%)")

full_PCA4
ggsave(filename = "full_PCA4.png", plot = full_PCA4, width = 8, height = 4, units = "in", dpi = 300)


