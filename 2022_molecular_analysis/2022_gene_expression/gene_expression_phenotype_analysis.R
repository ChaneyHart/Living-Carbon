#More thorough background analysis of comp between gene expression and physiology

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
library(PCAtools)
library(GGally)

#read in compiled csv
by_tree <- read.csv("Gene_exp_bytree_both_yrs.csv")
str(by_tree)
by_tree$PLGG1_2021 <- as.numeric(by_tree$PLGG1_2021)
by_tree$GDH_2021 <- as.numeric(by_tree$GDH_2021)
by_tree$MS_2021 <- as.numeric(by_tree$MS_2021)
by_tree$PLGG1_avg <- as.numeric(by_tree$PLGG1_avg)
by_tree$GDH_avg <- as.numeric(by_tree$GDH_avg)
by_tree$MS_avg <- as.numeric(by_tree$MS_avg)
by_tree$transgene_avg <- as.numeric(by_tree$transgene_avg)

#unfortunately only have data at the tree level in both years for a subset of the large block
#Trees in small block will not have any PLGG1 data from 2021 - omit them
by_tree_clean <- subset(by_tree, PLGG1_2021 > 0)
#For much of the analysis, only going to look at relationships amongst transgenic trees
transgene_by_tree <- subset(by_tree_clean, construct2 == "transgenic")



####################analysis option #1 ##########################

#run PCA on gene expression data inputting both years


PCA_opt1set <- transgene_by_tree[,6:11]
PCA_opt1 <- prcomp(PCA_opt1set, center = TRUE, scale. = TRUE)
summary(PCA_opt1)
#plot w/ vectors
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

#detach("package:reshape2", unload=TRUE)
#detach("package:reshape", unload=TRUE)

ggbiplot(PCA_opt1)
#cool, PLGG1 seems to be running in different direction than the two transgenes. But where are the events here?
#with IDs and event clusters
PCA_opt1_plot <- ggbiplot(PCA_opt1, obs.scale = 20, var.scale = 14, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (35.0%)") + ylab ("PC2 (32.1%)")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

PCA_opt1_plot
##
PCA_opt1_plot2 <- ggbiplot(PCA_opt1, obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short, var.axes=FALSE)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab ("PC1 (35.0%)") + ylab ("PC2 (32.1%)")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)+
  scale_color_discrete(labels=c('A','B','C','D','E','F','G','H'))

PCA_opt1_plot2
ggsave(filename = "PCA_opt1_plot2.png",plot = PCA_opt1_plot2, dpi=300)

##interpret PCA

PCA_opt1$rotation
#PLGG1 has contribution to PC1
#other transgenes contribute to PC2, especially GDH
cor(PCA_opt1set[,1:6])
# PLGG1 and GDH are indeed negatively correlated - though not very strongly
#Plot PLGG1 average against PC1



PCA_opt1_plot_34 <- ggbiplot(PCA_opt1, choices = c(3,4), obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

PCA_opt1_plot_34

#Which vars to include
plot(PCA_opt1)
screeplot(PCA_opt1, type = "line", main = "Scree plot")

PCA_opt1$sdev^2
# the first three are greater than 1. This is what I intuitively would have guessed as well

#Now need to assign the rotation values to individual trees and then summarize event means
opt1_PCA_score <- as.data.frame(PCA_opt1$x)

opt1_PCA_score$ID <- transgene_by_tree$ID
opt1_PCA_score$Event <- transgene_by_tree$Event

opt1_PCA_score_event <- opt1_PCA_score %>% group_by(Event)

opt1_gene_summary <- opt1_PCA_score_event %>% dplyr::summarise(
  n = n(),
  PC_Exp_1_sd = sd(PC1),
  PC_Exp_1 = mean(PC1),
  PC_Exp_1_se = (PC_Exp_1_sd/(sqrt(n))),
  PC_Exp_2_sd = sd(PC2),
  PC_Exp_2 = mean(PC2),
  PC_Exp_2_se = (PC_Exp_2_sd/(sqrt(n))),
  PC_Exp_3_sd = sd(PC3),
  PC_Exp_3 = mean(PC3),
  PC_Exp_3_se = (PC_Exp_3_sd/(sqrt(n)))
  
)


phys <- read.csv("../LC_June2022/Li6800_data/Aci_parameters_list.csv")
phys$ID <- phys$tree
physII <- read.csv("../LC_June2022/Li6800_data/Aci_parameters_list_flr.csv")
physII$ID <- physII$tree
growth <- read.csv("../LC_Nov_22_update/DBH_H_timeline_CT1_excluded_9_22.csv")
growth$VI <- growth$V419/1000

phys <- left_join(phys, growth, by = "ID")
physII <- left_join(physII, growth, by = "ID")

#find event means for phys and then growth variables
phys_event <- phys %>% group_by(event)
phys_summary <- phys_event %>% dplyr::summarise(
  n = n(),
  Ci_sd = sd(Ci_A_410),
  Ci = mean(Ci_A_410),
  Ci_se = (Ci_sd/(sqrt(n))),
  Amax_sd = sd(A_A_410),
  Amax = mean(A_A_410),
  Amax_se = (Amax_sd/(sqrt(n))),
  Vcmax_sd = sd(Vcmax),
  Vcmax = mean(Vcmax),
  Vcmax_se = (Vcmax_sd/(sqrt(n))),
  Jmax_sd = sd(Jmax),
  Jmax = mean(Jmax),
  Jmax_se = (Jmax_sd/(sqrt(n))),
  Jmax_Vcmax_ratio_sd = sd(Jmax_Vcmax_ratio),
  Jmax_Vcmax = mean(Jmax_Vcmax_ratio),
  Jmax_Vcmax_ratio_se = (Jmax_Vcmax_ratio_sd/(sqrt(n))),
  Rd_sd = sd(Rd),
  Rd = mean(Rd),
  Rd_se = (Rd_sd/(sqrt(n)))
)

phys_event_II <- physII %>% group_by(event)


physII_summary <- phys_event_II %>% dplyr::summarise(
  n = n(),
  PhiPS2_PhiCO2_sd = sd(PhiPS2_PhiCO2_A_410),
  PhiPS2_PhiCO2 = mean(PhiPS2_PhiCO2_A_410),
  PhiPS2_PhiCO2_se = (PhiPS2_PhiCO2_sd/(sqrt(n)))
)

phys_summary_III <- inner_join(phys_summary, physII_summary, by = "event")

growth_event <- growth %>% group_by(event)

growth_summary <- growth_event %>% dplyr::summarize(
  n = n(),
  VI_sd = sd(VI),
  VI = mean(VI),
  VI_se = (VI_sd/(sqrt(n))),
  Chl_sd = sd(Aug_SPAD),
  Chl = mean(Aug_SPAD),
  Chl_se = (Chl_sd/(sqrt(n)))
)

#join datasets of event means

growth_phys_summary <- inner_join(phys_summary_III, growth_summary, by = "event")

opt1_gene_summary$event <- opt1_gene_summary$Event
opt1_full_summary <- inner_join(opt1_gene_summary, growth_phys_summary, by = "event")

#
#examine relationship between variables via correlation matrix

opt1_full_corrset <- opt1_full_summary[,c(4,7,10,15,18,21,24,27,30,34,38,41)]
opt1_full_cor_plot <- ggpairs(opt1_full_corrset, upper = list(continuous = wrap("cor", size = 2)), lower = list(continuous = wrap("points", size = 0.5)), diag = list(continuous ="blankDiag"), axisLabels = "none")
ggsave(filename = "opt1_full_cor_plot.png", plot = opt1_full_cor_plot, width = 6, height = 4, dpi = 300)
#ok, pretty busy


opt1_full_set <- opt1_full_summary[,c(4,7,10,18,21,24,27,38,41)]
opt1_PCA_full <- prcomp(opt1_full_set, center = TRUE, scale. = TRUE)
summary(opt1_PCA_full)

ggbiplot(opt1_PCA_full)
#nice

opt1_full_summary$Event_short <- c("13-15B","13-15E","2H","4A","5A","5C","7")

opt1_full_PCAplot <- ggbiplot(opt1_PCA_full, obs.scale = 1.5, var.scale = 5, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("New_PC1 (40.5%)")+
  ylab("New_PC2 (28.4%)")+
  geom_point(aes(colour=opt1_full_summary$Event_short), size = 3)

opt1_full_PCAplot
ggsave(filename = "opt1_full_PCA.png", plot = opt1_full_PCAplot, width = 8, height = 4, units = "in", dpi = 300)

opt1_PCA_full$rotation
cor(opt1_full_set[,1:9])
plot(opt1_PCA_full)

opt1_full_PCAplot2 <- ggbiplot(opt1_PCA_full, choices = c(3,4),obs.scale = 1.5, var.scale = 1, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour=opt1_full_summary$Event_short), size = 3)
opt1_full_PCAplot2

#######
#plotting scatterplots to drive home interpretation of PCA
PC1_vol <- ggplot(opt1_full_summary, aes(x = PC_Exp_1, y = VI))+
  geom_point(aes(color = Event_short),size=4)+
  geom_errorbar(aes(xmin = PC_Exp_1-PC_Exp_1_se, xmax = PC_Exp_1+PC_Exp_1_se))+
  geom_errorbar(aes(ymin = VI-VI_se, ymax=VI+VI_se))+
  ylab("Mean stem volume index (cm^3)")
  
PC1_vol
ggsave(filename = "PC1_vol.png", plot = PC1_vol, dpi = 300)

PC_VI_model <- lm(VI~PC_Exp_1, opt1_full_summary)
summary(PC_VI_model)

PC2_vol <- ggplot(opt1_full_summary, aes(x = PC_Exp_2, y = VI))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin = PC_Exp_2-PC_Exp_2_se, xmax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(ymin = VI-VI_se, ymax=VI+VI_se))+
  ylab("Mean stem volume index (cm^3)")

PC2_vol
ggsave(filename = "PC2_vol.png", plot = PC2_vol, dpi = 300)
  

PC2_VI_model <- lm(VI~PC_Exp_2, opt1_full_summary)
summary(PC2_VI_model)


library(plotly)
plot_ly(opt1_full_summary, x = ~PC_Exp_1, y = ~PC_Exp_2, z = ~VI, color = ~Event_short, error_y = ~list(array = sd,color = '#000000'))
#add error bars

### another visualization option

ggplot(opt1_full_summary, aes(x=PC_Exp_1,y=PC_Exp_2,group = VI, color = VI))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(xmin = PC_Exp_1-PC_Exp_1_se, xmax = PC_Exp_1+PC_Exp_1_se))+
  geom_errorbar(aes(ymin = PC_Exp_2-PC_Exp_2_se, ymax = PC_Exp_2+PC_Exp_2_se))

###another
library(latticeExtra)

cloud(PC_Exp_1~PC_Exp_2+VI, opt1_full_summary, panel.3d.cloud=panel.3dbars, col.facet='grey', 
      xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
      par.settings = list(axis.line = list(col = "transparent")))

#####another####
install.packages("plot3D")
library("plot3D")

x <- opt1_full_summary$PC_Exp_1
y <- opt1_full_summary$PC_Exp_2
z <- opt1_full_summary$VI
#simple graph

scatter3D(x,y,z)

##create CI matrix
z_min <- (((opt1_full_summary$VI)-(opt1_full_summary$VI_se)))
z_max <- (((opt1_full_summary$VI)+(opt1_full_summary$VI_se)))
CI <- list(z = matrix(c(z_min,z_max), ncol =2))
CI$z

#build regression plane
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
#with regression surface
scatter3D(x,y,z,bty = "b2",pch = 19, cex = 2, ticktype = "detailed", col.var = z, theta = 25, phi =20, type = "h",clab = "Stem volume index (cm^3)",surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA, fit = fitpoints))

#simple
tiff('3d_plot.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
  scatter3D(x,y,z,bty = "b2",pch = 19, cex = 2, ticktype = "simple", col.var = z, theta =40, phi = 5, type = "h",clab = "Stem volume index (cm^3)")

dev.off()

text3D(x,y,z, labels = opt1_full_summary$Event_short,add = TRUE, colkey = TRUE, cex = 1)

threeD_p
#####Compare to photosynthetic efficiency########

ggplot(opt1_full_summary, aes(x = PC_Exp_1, y = Vcmax))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin = PC_Exp_1-PC_Exp_1_se, xmax = PC_Exp_1+PC_Exp_1_se))+
  geom_errorbar(aes(ymin = Vcmax-Vcmax_se, ymax=Vcmax+Vcmax_se))+
  ylab("mean Vcmax (µmol/m^s*s)")

ggplot(opt1_full_summary, aes(x = PC_Exp_1, y = Jmax))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin = PC_Exp_1-PC_Exp_1_se, xmax = PC_Exp_1+PC_Exp_1_se))+
  geom_errorbar(aes(ymin = Jmax-Jmax_se, ymax=Jmax+Jmax_se))+
  ylab("mean Jmax (µmol/m^s*s)")

PC1_Vcmax_model <- lm(Vcmax~PC_Exp_1, opt1_full_summary)
summary(PC1_Vcmax_model)

PC1_Jmax_model <- lm(Jmax~PC_Exp_1, opt1_full_summary)
summary(PC1_Jmax_model)


###PC2
ggplot(opt1_full_summary, aes(x = PC_Exp_2, y = Vcmax))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin = PC_Exp_2-PC_Exp_2_se, xmax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(ymin = Vcmax-Vcmax_se, ymax=Vcmax+Vcmax_se))+
  ylab("mean Vcmax (µmol/m^s*s)")

ggplot(opt1_full_summary, aes(x = PC_Exp_2, y = Jmax))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin = PC_Exp_2-PC_Exp_2_se, xmax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(ymin = Jmax-Jmax_se, ymax=Jmax+Jmax_se))+
  ylab("mean Jmax (µmol/m^s*s)")

PC2_Vcmax_model <- lm(Vcmax~PC_Exp_2, opt1_full_summary)
summary(PC2_Vcmax_model)

PC2_Jmax_model <- lm(Jmax~PC_Exp_2, opt1_full_summary)
summary(PC2_Jmax_model)

####################option 2 #########################

###### instead of putting gene expression PCs into full PCA, put raw gene expression data
transgene_by_tree

opt2_gene_by_event <- transgene_by_tree %>% group_by(Event)

opt2_gene_summary <- opt2_gene_by_event %>% dplyr::summarise(
  PLGG1_2021_sd = sd(PLGG1_2021),
  PLGG1_2021 = mean(PLGG1_2021),
  PLGG1_2021_n = n(), 
  PLGG1_2021_se = PLGG1_2021_sd / sqrt(PLGG1_2021_n),
  GDH_2021_sd = sd(GDH_2021),
  GDH_2021 = mean(GDH_2021),
  GDH_2021_n = n(), 
  GDH_2021_se = GDH_2021_sd / sqrt(GDH_2021_n),
  MS_2021_sd = sd(MS_2021),
  MS_2021 = mean(MS_2021),
  MS_2021_n = n(), 
  MS_2021_se = MS_2021_sd / sqrt(MS_2021_n),
  PLGG1_2022_sd = sd(PLGG1_2022),
  PLGG1_2022 = mean(PLGG1_2022),
  PLGG1_2022_n = n(), 
  PLGG1_2022_se = PLGG1_2022_sd / sqrt(PLGG1_2022_n),
  GDH_2022_sd = sd(GDH_2022),
  GDH_2022 = mean(GDH_2022),
  GDH_2022_n = n(), 
  GDH_2022_se = GDH_2022_sd / sqrt(GDH_2022_n),
  MS_2022_sd = sd(MS_2022),
  MS_2022 = mean(MS_2022),
  MS_2022_n = n(), 
  MS_2022_se = MS_2022_sd / sqrt(MS_2022_n),
  PLGG1_av_sd = sd(PLGG1_avg),
  PLGG1_av = mean(PLGG1_avg),
  PLGG1_av_se = (PLGG1_av_sd/(sqrt(PLGG1_2021_n))),
  MS_av_sd = sd(MS_avg),
  MS_n =n(),
  MS_av = mean(MS_avg),
  MS_av_se = MS_av_sd/(sqrt(MS_n)),
  GDH_av_sd = sd(GDH_avg),
  GDH_n =n(),
  GDH_av = mean(GDH_avg),
  GDH_av_se = GDH_av_sd/(sqrt(GDH_n)),
  transgene_av_sd = sd(transgene_avg),
  transgene_n = n(),
  transgene_av = mean(transgene_avg),
  transgene_av_se = transgene_av_sd/(sqrt(transgene_n))
)  

opt2_gene_summary$event <- opt2_gene_summary$Event
opt2_full_summary <- inner_join(opt2_gene_summary, growth_phys_summary, by = "event")

##

opt2_full_set <- opt2_full_summary[,c(3,7,11,15,19,23,32,35,38,41,44,52,55)]
opt2_PCA_full <- prcomp(opt2_full_set, center = TRUE, scale. = TRUE)
summary(opt2_PCA_full)

ggbiplot(opt2_PCA_full)
opt2_full_summary$Event_short <- c("13-15B","13-15E","2H","4A","5A","5C","7")


opt2_full_PCAplot <- ggbiplot(opt2_PCA_full, obs.scale = 2, var.scale = 3.5, ellipse = TRUE)+
  ggtitle("phenotype PCA by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  xlab("PC1 (40.5%)")+
  ylab("PC2 (28.4%")+
  geom_point(aes(colour=opt2_full_summary$Event_short), size = 3)

opt2_full_PCAplot

ggsave(filename = "opt2_full_PCAplot.png", plot = opt2_full_PCAplot, dpi = 300)

opt2_PCA_full$rotation
cor(opt2_full_set[,1:13])


### reality check looking at vars

opt2_full_summary$Event_short <- c("13-15B","13-15E","2H","4A","5A","5C","7")


VI_vcmax_event2 <- ggplot(opt1_full_summary,aes(y = VI, x =Vcmax))+
  geom_point(aes( color = Event_short), size = 3)+
  ylab('volume index (mm3)')+
  xlab('Vcmax (µmol/mol s)')+
  geom_errorbar(aes(xmin = Vcmax-Vcmax_se, xmax = Vcmax+Vcmax_se))+
  geom_errorbar(aes(ymin = VI-VI_se, ymax = VI+VI_se))
VI_vcmax_event2
ggsave(filename = "VI_vcmax_event2.png", plot = VI_vcmax_event2, dpi = 300)

VI_model2 <- lm(VI~Vcmax, opt1_full_summary)
summary(VI_model)




VI_Jmax_event2 <- ggplot(opt1_full_summary,aes(y = VI, x =Jmax))+
  geom_point(aes(color = Event_short), size = 3)+
  ylab('volume index (mm3)')+
  xlab('Jmax (µmol/mol s)')+
  geom_errorbar(aes(xmin = Jmax-Jmax_se, xmax = Jmax+Jmax_se))+
  geom_errorbar(aes(ymin = VI-VI_se, ymax = VI+VI_se))
VI_Jmax_event2
ggsave(filename = "VI_Jmax_event2.png", plot = VI_Jmax_event2, dpi = 300)

VI_model2 <- lm(VI~Jmax, opt1_full_summary)
summary(VI_model2)

growth_phys_model <- lm(V419~Vcmax, Aci_growth)
growth_phys_modelII <- lm(V419~Jmax, Aci_growth)
summary(growth_phys_model)
summary(growth_phys_modelII)
# r2 near 0
growth_phys


cor(growth_phys_summary[,c(4,7,10)])

###
#plot PLGG1 against volume
ggplot(opt2_full_summary, aes(x=PLGG1_av, y = VI))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin=PLGG1_av-PLGG1_av_se, xmax = PLGG1_av+PLGG1_av_se))+
  geom_errorbar(aes(ymin=VI-VI_se, ymax=VI+VI_se))+
  xlab("Average PLGG1 expression (rel. fold change)")+
  ylab("stem volume index (mm^3)")

PLGG1_VI_model <- lm(VI~PLGG1_av, opt2_full_summary)
summary(PLGG1_VI_model)


###plot transgenes against volume
ggplot(opt2_full_summary, aes(x=GDH_2022, y = VI))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin=PLGG1_av-PLGG1_av_se, xmax = PLGG1_av+PLGG1_av_se))+
  geom_errorbar(aes(ymin=VI-VI_se, ymax=VI+VI_se))+
  xlab("Average GDH expression (rel. fold change)")+
  ylab("stem volume index (mm^3)")
####
#include controls in plotting PLGG1 with volume. Will have to re-summarize by events

opt2_2_gene_by_event <- by_tree %>% group_by(Event)

opt2_2_gene_summary <- opt2_2_gene_by_event %>% dplyr::summarise(
  PLGG1_2021_sd = sd(PLGG1_2021),
  PLGG1_2021 = mean(PLGG1_2021),
  PLGG1_2021_n = n(), 
  PLGG1_2021_se = PLGG1_2021_sd / sqrt(PLGG1_2021_n),
  GDH_2021_sd = sd(GDH_2021),
  GDH_2021 = mean(GDH_2021),
  GDH_2021_n = n(), 
  GDH_2021_se = GDH_2021_sd / sqrt(GDH_2021_n),
  MS_2021_sd = sd(MS_2021),
  MS_2021 = mean(MS_2021),
  MS_2021_n = n(), 
  MS_2021_se = MS_2021_sd / sqrt(MS_2021_n),
  PLGG1_2022_sd = sd(PLGG1_2022),
  PLGG1_2022 = mean(PLGG1_2022),
  PLGG1_2022_n = n(), 
  PLGG1_2022_se = PLGG1_2022_sd / sqrt(PLGG1_2022_n),
  GDH_2022_sd = sd(GDH_2022),
  GDH_2022 = mean(GDH_2022),
  GDH_2022_n = n(), 
  GDH_2022_se = GDH_2022_sd / sqrt(GDH_2022_n),
  MS_2022_sd = sd(MS_2022),
  MS_2022 = mean(MS_2022),
  MS_2022_n = n(), 
  MS_2022_se = MS_2022_sd / sqrt(MS_2022_n),
  PLGG1_av_sd = sd(PLGG1_avg),
  PLGG1_av = mean(PLGG1_avg),
  PLGG1_av_se = (PLGG1_av_sd/(sqrt(PLGG1_2021_n)))
  
)  

opt2_2_gene_summary$event <- opt2_2_gene_summary$Event
opt2_2_full_summary <- inner_join(opt2_2_gene_summary, growth_phys_summary, by = "event")

##plot
ggplot(opt2_2_full_summary, aes(x=PLGG1_av, y = VI))+
  geom_point(aes(color = Event_short),size =4)+
  geom_errorbar(aes(xmin=PLGG1_av-PLGG1_av_se, xmax = PLGG1_av+PLGG1_av_se))+
  geom_errorbar(aes(ymin=VI-VI_se, ymax=VI+VI_se))+
  xlab("Average PLGG1 expression (rel. fold change)")+
  ylab("stem volume index (mm^3)")

PLGG1_VI_model_2 <- lm(VI~PLGG1_av, opt2_2_full_summary)
summary(PLGG1_VI_model)


###PCA option 3
###plot PLGG1 2021 and PLGG1 2022 separately

PCA_opt3set1 <- transgene_by_tree[,c(6,7,8,10,11)]
PCA_opt3set2 <- transgene_by_tree[,c(7,8,9,10,11)]
###

PCA_opt3_1 <- prcomp(PCA_opt3set1, center = TRUE, scale. = TRUE)
PCA_opt3_2 <- prcomp(PCA_opt3set2, center = TRUE, scale. = TRUE)
summary(PCA_opt3_1)
summary(PCA_opt3_2)
#plot w/ vectors
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

#detach("package:reshape2", unload=TRUE)
#detach("package:reshape", unload=TRUE)

ggbiplot(PCA_opt3_1)
#cool, PLGG1 seems to be running in different direction than the two transgenes. But where are the events here?
#with IDs and event clusters
PCA_opt3_1_plot <- ggbiplot(PCA_opt3_1, obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

PCA_opt3_1_plot
###
ggbiplot(PCA_opt3_2)
#cool, PLGG1 seems to be running in different direction than the two transgenes. But where are the events here?
#with IDs and event clusters
PCA_opt3_2_plot <- ggbiplot(PCA_opt3_2, obs.scale = 1, var.scale = 1, labels.size = 2, varname.size = 3, ellipse = TRUE, groups = transgene_by_tree$Event_short)+
  ggtitle("PCA of field gene expression by event")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  geom_point(aes(colour= transgene_by_tree$Event_short), size = 1)

PCA_opt3_2_plot


###plotting genes against PCs for PAG presentation
full_comp <- inner_join(opt1_full_summary, opt2_full_summary, by = "Event")
PC1_PLGG1 <- ggplot(full_comp, aes(x=PLGG1_av, y = PC_Exp_1))+
  geom_point(aes(color = Event_short.x, size =4))+
  geom_errorbar(aes(ymin = PC_Exp_1-PC_Exp_1_se, ymax = PC_Exp_1+PC_Exp_1_se))+
  geom_errorbar(aes(xmin = PLGG1_av-PLGG1_av_se, xmax = PLGG1_av+PLGG1_av_se))
PC1_PLGG1

ggsave(filename = "PC1_PLGG1.png", plot = PC1_PLGG1, dpi = 300)

ggplot(full_comp, aes(x=GDH_av, y = PC_Exp_2))+
  geom_point(aes(color = Event_short.x, size =4))+
  geom_errorbar(aes(ymin = PC_Exp_2-PC_Exp_2_se, ymax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(xmin = GDH_av-GDH_av_se, xmax = GDH_av+GDH_av_se))

ggplot(full_comp, aes(x=MS_av, y = PC_Exp_2))+
  geom_point(aes(color = Event_short.x, size =4))+
  geom_errorbar(aes(ymin = PC_Exp_2-PC_Exp_2_se, ymax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(xmin = MS_av-MS_av_se, xmax = MS_av+MS_av_se))

transgene_PC2 <- ggplot(full_comp, aes(x=transgene_av, y = PC_Exp_2))+
  geom_point(aes(color = Event_short.x, size =4))+
  geom_errorbar(aes(ymin = PC_Exp_2-PC_Exp_2_se, ymax = PC_Exp_2+PC_Exp_2_se))+
  geom_errorbar(aes(xmin = transgene_av-transgene_av_se, xmax = transgene_av+transgene_av_se))

ggsave(filename = "transgene_PC2.png", plot = transgene_PC2, dpi = 300)
####
#creating graph for PAG presentation
transgene_by_tree_av <- transgene_by_tree[,c(1,3,12,13,14)]
#install.packages("reshape2")                                 
library("reshape2") 
transgene_by_tree_av_long <- melt(transgene_by_tree_av)
head(transgene_by_tree_av_long)
gene_raw <- ggplot(transgene_by_tree_av_long, aes(x = variable, y =value, fill = variable))+
  geom_boxplot()+
  geom_point(aes(x= variable, y=value, color = Event_short))

gene_raw
ggsave(filename = "gene_raw.png", plot = gene_raw, width = 7, height = 3, units = "in", dpi = 300)
