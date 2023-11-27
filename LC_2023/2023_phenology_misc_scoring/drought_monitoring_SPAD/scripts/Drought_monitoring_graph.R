#Drought_monitoring_graph




#########
## Spatial changes in drought scores
#read in  cleaned data
drought_scores_2023 <- read.csv(file = "LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/Drought_scores_clean_2023.csv")
#install.packages("latticeExtra")
library(latticeExtra)

#drought effect varied by position
drought_day1 <- levelplot(DS_1 ~-column*row, drought_scores_2023,panel = panel.levelplot.points, cex =1.2, row.values = drought_scores_2023$row, column.values = drought_scores_2023$column)+
  layer_(panel.2dsmoother(..., n=424))
drought_day1
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_day1_spatial.png", plot = drought_day1,dpi=300)

drought_day18 <- levelplot(DS_18 ~-column*row, drought_scores_2023,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
drought_day18
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_day18_spatial.png", plot = drought_day18,dpi=300)


drought_day32 <- levelplot(DS_32 ~-column*row, drought_scores_2023,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
drought_day32
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_day32_spatial.png", plot = drought_day32,dpi=300)


drought_day55 <- levelplot(DS_55 ~-column*row, drought_scores_2023,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
drought_day55
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_day55_spatial.png", plot = drought_day55,dpi=300)

drought_scores_2023$delta_drought <- drought_scores_2023$DS_55 - drought_scores_2023$DS_1

delta_drought <- levelplot(delta_drought ~-column*row, drought_scores_2023,panel = panel.levelplot.points, cex =1.2)+
  layer_(panel.2dsmoother(..., n=424))
delta_drought
ggsave(filename = "LC_2023/2023_phenology_misc_scoring/drought_delta_spatial.png", plot = drought_day55,dpi=300)


#How did drought scores change over time

hist(drought_scores_2023$DS_1,xlab = "drought score")
hist(drought_scores_2023$DS_18,xlab = "drought score")
hist(drought_scores_2023$DS_32,xlab = "drought score")
hist(drought_scores_2023$DS_55,xlab = "drought score")
str(drought_scores_2023)
drought_scores_2023$DS_1 <- as.numeric(drought_scores_2023$DS_1)


##############

drought_summary_long <- read.csv("LC_2023/2023_phenology_misc_scoring/drought_monitoring_SPAD/clean/drought_summary_long.csv")

myColors <- c('gray',
              'cyan4',
              'red4')

names(myColors) <- levels(drought_summary_long$construct)

#graph day 1
drought_summary_by_day <- drought_summary_long %>% group_by(Day) %>% summarize(
  percent_0 = mean(percent_0),
  percent_1 = mean(percent_1),
  percent_2 = mean(percent_2),
  percent_3 = mean(percent_3),
  percent_4 = mean(percent_4),
  percent_5 = mean(percent_5)
)

ggplot(subset(drought_summary_by_day), aes(x=Day))+
  geom_point(aes(y=percent_0))+
  geom_point(aes(y=percent_1))+
  geom_point(aes(y=percent_3))+
  geom_point(aes(y=percent_4))+
  geom_line(aes(y=percent_0),color="red3")+
  geom_line(aes(y=percent_1),color="blue3")+
  geom_line(aes(y=percent_3),color="purple3")+
  geom_line(aes(y=percent_4),color="orange")+
  ylab("Fraction of trees")+
  xlab("Days since June 20th")
  
percent_0 <- ggplot(drought_summary_long, aes(x=Day,y=percent_0,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of non-drought trees (recieved score 0)")+
  theme_bw()

percent_0 + colorScale

percent_1 <- ggplot(drought_summary_long, aes(x=Day,y=percent_1,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of yellowing trees (recieved score 1)")+
  theme_bw()

percent_1 + colorScale

percent_2 <- ggplot(drought_summary_long, aes(x=Day,y=percent_2,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of trees with leaf drop (recieved score 2)")+
  theme_bw()

percent_2 + colorScale


percent_3 <- ggplot(drought_summary_long, aes(x=Day,y=percent_3,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of trees with substantial leaf drop (recieved score 3)")+
  theme_bw()

percent_3 + colorScale

percent_4 <- ggplot(drought_summary_long, aes(x=Day,y=percent_4,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of trees with complete leaf drop (recieved score 4)")+
  theme_bw()

percent_4 + colorScale


percent_5 <- ggplot(drought_summary_long, aes(x=Day,y=percent_5,group=event,color=construct, shape = event))+
  geom_line(aes(color=construct))+
  geom_point(aes(shape=event))+
  xlab("Days since July 20th 2023")+
  ylab("proportion of dead trees (recieved score 5)")+
  theme_bw()

percent_5 + colorScale


