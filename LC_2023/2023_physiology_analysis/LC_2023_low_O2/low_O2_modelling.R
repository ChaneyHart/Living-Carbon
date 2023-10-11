#low O2 A-ETR modelling


library(dplyr)
library(ggplot2)

low_O2_dat_1 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_low_O2/low_O2_1.csv")
low_O2_dat_2 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_low_O2/low_O2_2.csv")
low_O2_dat_2 <- subset(low_O2_dat_2, ID != "LCOR-310")
low_O2_dat_3 <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_low_O2/low_O2_3.csv")

low_O2_dat <- rbind(low_O2_dat_1,low_O2_dat_2,low_O2_dat_3)
low_O2_dat$PhiPS2 <- as.numeric(low_O2_dat$PhiPS2)
low_O2_dat$ETR <- as.numeric(low_O2_dat$ETR)


low_O2_dat <- subset(low_O2_dat, A > 0)


low_02_plot <- ggplot(low_O2_dat, aes(x=PhiCO2, y=PhiPS2, color = Event_short))+
  geom_point()

ggsave(filename = "LC_2023/2023_physiology_analysis/LC_2023_low_O2/low_O2_plot.png", plot = low_02_plot, dpi = 300)

Low_O2_fit <- lm(PhiPS2~PhiCO2, low_O2_dat)
summary(Low_O2_fit)

low_O2_fit_yint <- Low_O2_fit$coefficients[1]
low_O2_fit_coeff1 <- Low_O2_fit$coefficients[2]

