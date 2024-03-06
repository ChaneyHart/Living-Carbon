#gene expression PCA, 2024

library(dplyr)
library(ggplot2)
library(gridExtra)
library((ggthemes))
library(ggsignif)
library(GGally)
library(ggbiplot)
library(tidyverse)

gene_exp <- read.csv(file = "LC_2023/2023_molecular_analysis/gene_exp.csv")

gene_exp_summary <- gene_exp %>% group_by(Event,sample_date,gene) %>% summarize(
  n = n(),
  std = sd(expression_level),
  mean = mean(expression_level),
  se = std/(sqrt(n))
)

gene_exp_summary <- gene_exp_summary %>% mutate(
  Class = case_when(
    Event == "5A" | Event == "5C" | Event == "4A" ~ "intermediate",
    Event == "13-15E" | Event == "2H" ~ "high",
    Event == "16-20" | Event == "8-9D" ~ "Control",
    Event == "CT3" ~ "WT"))


gene_exp_summary_wide <- pivot_wider(gene_exp_summary, id_cols = c(Event,Class),names_from = c("gene","sample_date"),values_from = c("mean","se"))

#plot to explore
gene_exp_summary$sample_date <- as.factor(gene_exp_summary$sample_date)
gene_exp_summary$gene <- as.factor(gene_exp_summary$gene)

mycolors <- c("gray1","chartreuse4","indianred3","gray")

gene_exp_summary <- gene_exp_summary %>% mutate(
  sample_date = factor(sample_date,
                       levels = c("09-12-2021","06-20-2022","07-14-2023","08-14-2023")
))
gene_expression_yearly <- ggplot(gene_exp_summary, aes(x=gene, y=mean, color=Class))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.1,alpha=0.5)+
  scale_color_manual(values = mycolors)+
  theme_bw()+
  ylab("mean relative expression")+
  facet_grid(~sample_date)

ggsave(filename="LC_2023/2023_molecular_analysis/gene_expression_yearly.png",plot=gene_expression_yearly, dpi=300)

ggplot(gene_exp_summary, aes(x=gene,y=mean,color=Class, shape=sample_date))+
  geom_jitter(width=0.2)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.1,alpha=0.5)+
  scale_color_manual(values = mycolors)+
  theme_bw()


PCA_set <- gene_exp_summary_wide[c(1,3,4,5,6),3:14]
colnames(PCA_set) <- c("GDH_22","MS_22","PLGG1_22","GDH_7_23","MS_7_23","PLGG1_7_23","GDH_8_23","MS_8_23","PLGG1_8_23","GDH_21","MS_21","PLGG1_21")
PCA_event <- prcomp(PCA_set, center = TRUE, scale. = TRUE)
summary(PCA_event)
#plot w/ vectors
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

#detach("package:reshape2", unload=TRUE)
#detach("package:reshape", unload=TRUE)

PCA_colors <- c("")
ggbiplot(PCA_event)
mycolors_pca <- c("chartreuse4","indianred3")
#cool, PLGG1 seems to be running in different direction than the two transgenes. But where are the events here?
#with IDs and event clusters
event_PCA <- ggbiplot(PCA_event, obs.scale = 1, var.scale = 3,labels.size = 0.2)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (68.9%)") + ylab ("PC2 (18.5%)")+
  geom_point(aes(colour=c("high","high","intermediate","intermediate","intermediate"), size = 1))+
  scale_color_manual(values = mycolors_pca)
  


event_PCA
ggsave(filename = "LC_2023/2023_molecular_analysis/event_PCA.png", plot = event_PCA, width = 12, height = 4, units = "in",dpi = 300)

event_PCA2 <- ggbiplot(PCA_event, choice = c(2,3),obs.scale = 1, var.scale = 1,labels.size = 0.2)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=c("high","high","intermediate","intermediate","intermediate"), size = 1))+
  scale_color_manual(values = mycolors_pca)



event_PCA2
ggsave(filename = "LC_2023/2023_molecular_analysis/event_PCA2.png", plot = event_PCA2, width = 12, height = 4, units = "in",dpi = 300)



event_PCA3 <- ggbiplot(PCA_event, choice = c(1,3),obs.scale = 1, var.scale = 1,labels.size = 0.2)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=c("high","high","intermediate","intermediate","intermediate"), size = 1))+
  scale_color_manual(values = mycolors_pca)



event_PCA3
ggsave(filename = "LC_2023/2023_molecular_analysis/event_PCA3.png", plot = event_PCA3, width = 12, height = 4, units = "in",dpi = 300)
