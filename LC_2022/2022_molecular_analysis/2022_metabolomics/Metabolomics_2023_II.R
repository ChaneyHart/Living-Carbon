##Metabolics 2023_revised_6_7_23
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(dplyr)
library(ggsignif)
library(PCAtools)
library(GGally)
library(ggbiplot)



metabolomics <- read.csv(file = "2022_molecular_analysis/2022_metabolomics/Meabolomics_2023_normalized.csv")
str(metabolomics)
metabolomics$file_name
metabolomics$ID <- substr(metabolomics$file_name,1,8)
#remove quality check samples
metabolomics <- metabolomics[-1:-16,]
#remove Ribitol, the internal standard
metabolomics <- metabolomics[,-2]
#remove columns with NA
metabolomics <- metabolomics[ , colSums(is.na(metabolomics)) == 0]


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
  geom_signif(comparisons = list(c("2H","CT3")),map_signif_level=TRUE, y_position = 0.001)+
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


ggplot(metabolomics_IV, aes(x=V497, y=glycolic.acid))+
  geom_point(aes(colour = event_short))

######

epichatechin_event <- ggplot(metabolomics_III, aes(x=event_short, y = X....epicatechin, fill = construct2.x))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("13-15E","4A")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("2H","4A")),map_signif_level=TRUE, y_position = 0.001)+
  ylab("Relative epicatechin level (to IS)")+
  xlab("event")

epichatechin_event

chlorogenic_acid_event <- ggplot(metabolomics_III, aes(x=event_short, y = chlorogenic.acid.2, fill = construct2.x))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("13-15E","4A")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("2H","4A")),map_signif_level=TRUE, y_position = 0.001)+
  ylab("Relative chlorogenic acid level (to IS)")+
  xlab("event")

chlorogenic_acid_event


####PCA of metabolites ########

##part of purpose is to explore feature reduction

#subset to get metabolites and PLGG1_2022

PCA_metabolites_set <- metabolomics_III[,c(2:61,67,111)]

PCA_meta <- subset(metabolomics_III, select = c("ID","event_short","block","PLGG1_2022","GDH_2022","MS_2022"))
#exporting CsV for python
PCA_metabolites_export <- inner_join(PCA_metabolites_set, PCA_meta, by = "ID")
write.csv(PCA_metabolites_export, file = "PCA_metabolites_export.csv")
##


PCA_metabolites <- prcomp(PCA_metabolites_set[,c(1:59)], center = TRUE, scale. = TRUE, print = TRUE)


#which PCs are most important?
ggscreeplot(PCA_metabolites)

#looks like first 11 are useful

ggbiplot(PCA_metabolites)
#The main distinguisher is that some trees seem to have less of everything
#sucrose is set off to the side as well
ggbiplot(PCA_metabolites, choices = c(1,3))
ggbiplot(PCA_metabolites, choices = c(1,4))
ggbiplot(PCA_metabolites, choices = c(2,3))
ggbiplot(PCA_metabolites, choices = c(2,4))

##how do different events compare in PC1 versus PC2 space
metabolite_PCA <- ggbiplot(PCA_metabolites,obs.scale = 1, var.axes = FALSE, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_metabolites_set$event_short )+
  ggtitle("Principal components of poplar leaf metabolome")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (14.7%)") + ylab ("PC2 (11.9%)")+
  geom_point(aes(colour= PCA_metabolites_set$event_short), size = 1)

ggsave(filename = "2022_molecular_analysis/2022_metabolomics/metabolite_PCA.png",plot = metabolite_PCA, dpi = 300)
#no one event is clustered away from the rest. 

PCA_2_3 <- ggbiplot(PCA_metabolites,obs.scale = 1, choices = c(2,3), var.axes = FALSE, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_metabolites_set$event_short )+
  ggtitle("Principal components of poplar leaf metabolome")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC2 (11.9%)") + ylab("PC3 (9.7%)")+
  geom_point(aes(colour= PCA_metabolites_set$event_short), size = 1)


PCA_2_3
ggsave(filename = "2022_molecular_analysis/2022_metabolomics/metabolite_PCA_2_3.png",plot = PCA_2_3, dpi = 300)

PCA_3_4 <- ggbiplot(PCA_metabolites,obs.scale = 1, choices = c(3,4), var.axes = FALSE, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = PCA_metabolites_set$event_short )+
  ggtitle("Principal components of poplar leaf metabolome")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC3 (9.7%)") + ylab("PC4 (6.9%)")+
  geom_point(aes(colour= PCA_metabolites_set$event_short), size = 1)

PCA_3_4
ggsave(filename = "2022_molecular_analysis/2022_metabolomics/metabolite_PCA_3_4.png", plot = PCA_3_4, dpi = 300)
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

PC1_min_graph <- ggplot(subset(PC1_loadings, PC1 < 0.05), aes(x=reorder(metabolite, PC1), y = PC1))+
  geom_bar(stat='identity', fill = "orange3")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC1 loading")+
  xlab("metabolite")

PC1_min_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC1_loadings_min.png",plot = PC1_min_graph,  height = 5, width = 4, units = "in", dpi = 300)

PC1_max_graph <- ggplot(subset(PC1_loadings, PC1 > 0.05), aes(x=reorder(metabolite, -1*PC1), y = PC1))+
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

####
#PC3
PC4_max <- slice_max(PCA_loadings,order_by = PC4, n= 10)
PC4_max <- subset(PC4_max, select = c("metabolite","PC4"))
PC4_min <- slice_min(PCA_loadings,order_by = PC4, n=10)
PC4_min <- subset(PC4_min, select = c("metabolite","PC4"))

PC4_loadings <- rbind(PC4_max, PC4_min)

PC4_min_graph <- ggplot(subset(PC4_loadings, PC4 < 0), aes(x=reorder(metabolite, PC4), y = PC4))+
  geom_bar(stat='identity', fill = "orange4")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC4 loading")+
  xlab("metabolite")

PC4_min_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC4_loadings_min.png",plot = PC4_min_graph,  height = 5, width = 4, units = "in", dpi = 300)


PC4_max_graph <- ggplot(subset(PC4_loadings, PC4 > 0), aes(x=reorder(metabolite, -1*PC4), y = PC4))+
  geom_bar(stat='identity', fill = "blue4")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("PC4 loading")+
  xlab("metabolite")


PC4_max_graph 
ggsave(file = "2022_molecular_analysis/2022_metabolomics/PC4_loadings_max.png",plot = PC4_max_graph, height = 5, width = 4, units = "in", dpi = 300)




##assign event level mean and se of PC1-11
metabolites_PCA_score <- as.data.frame(PCA_metabolites$x)

metabolites_PCA_score$ID <- PCA_metabolites_set$ID
#Export individual tree level PCA stats for each tree
write.csv(metabolites_PCA_score,file="2022_molecular_analysis/2022_metabolomics/metabolome_PCA_scores.csv")

metabolites_PCA_score$event_short <- PCA_metabolites_set$event_short

metabolites_PCA_score_event <- metabolites_PCA_score %>% group_by(event_short)

metabolites_event_summary <- metabolites_PCA_score_event %>% dplyr::summarise(
  n = n(),
  PC_Exp_1_sd = sd(PC1),
  PC_Exp_1 = mean(PC1),
  PC_Exp_1_se = (PC_Exp_1_sd/(sqrt(n))),
  PC_Exp_2_sd = sd(PC2),
  PC_Exp_2 = mean(PC2),
  PC_Exp_2_se = (PC_Exp_2_sd/(sqrt(n))),
  PC_Exp_3_sd = sd(PC3),
  PC_Exp_3 = mean(PC3),
  PC_Exp_3_se = (PC_Exp_3_sd/(sqrt(n))),
  PC_Exp_4_sd = sd(PC4),
  PC_Exp_4 = mean(PC4),
  PC_Exp_4_se = (PC_Exp_4_sd/(sqrt(n))),
  PC_Exp_5_sd = sd(PC5),
  PC_Exp_5 = mean(PC5),
  PC_Exp_5_se = (PC_Exp_5_sd/(sqrt(n))),
  PC_Exp_6_sd = sd(PC6),
  PC_Exp_6 = mean(PC6),
  PC_Exp_6_se = (PC_Exp_6_sd/(sqrt(n))),
  PC_Exp_7_sd = sd(PC7),
  PC_Exp_7 = mean(PC7),
  PC_Exp_7_se = (PC_Exp_7_sd/(sqrt(n))),
  PC_Exp_8_sd = sd(PC8),
  PC_Exp_8 = mean(PC8),
  PC_Exp_8_se = (PC_Exp_8_sd/(sqrt(n))),
  PC_Exp_9_sd = sd(PC9),
  PC_Exp_9 = mean(PC9),
  PC_Exp_9_se = (PC_Exp_9_sd/(sqrt(n))),
  PC_Exp_10_sd = sd(PC10),
  PC_Exp_10 = mean(PC10),
  PC_Exp_10_se = (PC_Exp_10_sd/(sqrt(n))),
  PC_Exp_11_sd = sd(PC11),
  PC_Exp_11 = mean(PC11),
  PC_Exp_11_se = (PC_Exp_11_sd/(sqrt(n))))







###
#Classify trees based on PLGG1 expression levels

for (x in fruits) {
  if (x == "cherry") {
    break
  }
  print(x)
}
PCA_metabolites_set

PCA_metabolites_set$Class <- numeric(118)                       # Create empty data object
PCA_metabolites_set$Class[PCA_metabolites_set$PLGG1_2022 < 0.4] <- 1                   # Assign categories based on numeric range
PCA_metabolites_set$Class[PCA_metabolites_set$PLGG1_2022 > 0.8] <- 2

Metabolites_PLGG1_classes <- subset(PCA_metabolites_set, Class > 0)

Metabolite_PLGG1_export <- Metabolites_PLGG1_classes[,c(-62,-61,-60)]

Metabolite_PLGG1_export <- Metabolite_PLGG1_export %>% select(Class, everything())
#Metabolite_PLGG1_export <- Metabolite_PLGG1_export %>% select(ID, everything())

write.csv(Metabolite_PLGG1_export, file = "Metabolite_PLGG1_export.csv")


#down sample for equal class distributions
count(Metabolite_PLGG1_export$Class)

Metabolite_PLGG1_export <- Metabolite_PLGG1_export %>% group_by(Class) %>% sample_n(35)



# fit lm for each metabolite
storage <- list()
model <- list()
sig <- list()

for(i in names(Metabolite_PLGG1_export)[-1]){
  storage[[i]] <- lm(Class ~ get(i), Metabolite_PLGG1_export)
  sig[[i]] <- summary(lm(Class ~ get(i), Metabolite_PLGG1_export))$coefficients[,4][2]
}


sigs <- as.data.frame(sig[[,]][2])


sig[[1]][1]



summary(fit)$coefficients[,4]  






for(i in names(Metabolite_PLGG1_export)[-1]){
  storage[[i]] <- lm(get(i) ~ Class, Metabolite_PLGG1_export)
}

for(i in names(Metabolite_PLGG1_export[,-1])){
  model = lm(i~Class, data=Metabolite_PLGG1_export)
}


storage[[6]][1]




lm(Class ~ glycolic.acid, data = Metabolite_PLGG1_export)
summary(lm(Class ~ glycolic.acid, data = Metabolite_PLGG1_export))

##Read in photosynthetic metadata

dat_6_24 <- read.csv(file = "../LC_June2022/June_RNA_metabolomics/li600_raw_RNA_6_24/2022-06-24/RNA_li600_6_24_raw.csv",skip=1)



#Calculate sample size needed to determine difference in glycolate levels

#What is standard deviation in glycolate levels
mean_glycolate <- mean(metabolomics_III$glycolic.acid)
sd_glycolate <- sd(metabolomics_III$glycolic.acid)
var_glycolate <- var(metabolomics_III$glycolic.acid)
qt(0.025,27)
CV_glycolate <- mean_glycolate/sd_glycolate

Tval <- -1*qt(0.025,42)

#We could estimate a difference of this value as significant with 8 samples
detectable_difference <- Tval*(sqrt(((2*sd_glycolate^2)/8)))

((Tval^2)*(1.23^2))/(0.8^2)




#last year the difference between 2H and CT3 and 13-15E was:
metabolomics_event <- metabolomics_III %>% group_by(event_short)

metabolites_glycolate_event <- metabolomics_event %>% dplyr::summarise(
  glycolate = mean(glycolic.acid))

#2H and CT3

metabolites_glycolate_event[5,2] - metabolites_glycolate_event[11,2]

#13-15E and CT3

metabolites_glycolate_event[2,2] - metabolites_glycolate_event[11,2]


#To detect a difference of 0.00050000, how many samples do we need?

sqrt.n = (((2*Tval)*(sd_glycolate^2)))/0.0005000
n <- sqrt.n^2



sqrt_n_2 <- ((Tval*var_glycolate)/(0.00001)
n_2 <- sqrt_n_2^2

var_glycolate/0.0001
