#Process A-Ci data

#read in compiled datasets

ACI_compiled <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/CO2_response_curves_compiled.csv")
ACI_compiled$PhiPS2 <- as.numeric(ACI_compiled$PhiPS2)
ACI_compiled$ETR <- as.numeric(ACI_compiled$ETR)

#identifying outliers
hist(ACI_compiled$A)
hist(ACI_compiled$Ci)
qplot(Ci,A, data = ACI_compiled)
#some outliers to corrct
outlier_check <- lm(A~Ci, data = ACI_compiled)
plot(outlier_check, which = 1)

#first outlier identified point 420, LCOR 204
ACI_compiled[420,]
#imput assumed A based on model
#fx for viewing data
viewdat = function(dat){plot(y = dat$A, x = dat$Ci);abline(h = seq(0,20,1));dat}

LCOR204 = subset(ACI_compiled, ID == "LCOR-204")
viewdat(LCOR204)
LCOR204 <- LCOR204[-15,]
viewdat(LCOR204)
LCOR204_fit = fitaci(LCOR204,Tcorrect = "FALSE",varnames = list(ALEAF = 'A',Tleaf = 'Tleaf',Ci = 'Ci',PPFD = 'Qin'))
LCOR204_fit$Photosyn(Ci=2000)
#use fit to inform A
ACI_compiled[420,9] <- 26.09878
ACI_compiled[420,11] <- 2000*0.6667

#next outliers, Ci is negative. LCOR 502 (obs 547), LCOR-577 (obs 246), LCOR-164 (obs 610) LCOR-516 (obs(486))

ACI_neg_ci <- subset(ACI_compiled, Ci < 0)
neg_ci_names <- row.names(ACI_neg_ci)

#Imput Ci from Ca times 0.6667
for (i in neg_ci_names) {
  print(i)
  ACI_compiled[neg_ci_names, 11] <- ACI_compiled[neg_ci_names, 10]*0.666
}

#repeat outlier check

outlier_check <- lm(A~Ci, data = ACI_compiled)
plot(outlier_check, which = 1)
#Assumption of equal variance looks okaaaay, not great
plot(outlier_check, which = 2)
#normality not bad
plot(outlier_check, which = 3)
#nothing with standardized residual beyond 2 or -2
plot(outlier_check, which = 4)
plot(outlier_check, which = 5)

#another Ci to fix
ACI_compiled[280,11] <- 400*0.667



#need to first average assimilation, gsw and Ci for cases where there are multiple measurements at same CO2 setpoint
ACI_summary <- ACI_compiled %>% group_by(ID,Csetpoint) %>% summarise(
  E = mean(E),
  A = mean(A),
  Ci = mean(Ci),
  Qin = mean(Qin),
  gsw = mean(gsw),
  VPDlead = mean(VPDleaf),
  Fs = mean(Fs),
  Fm_prime = mean(Fm.),
  PhiPS2 = mean(PhiPS2),
  PhiCO2 = mean(PhiCO2),
  ETR = mean(ETR),
  Tleaf = mean(Tleaf),
  sample_date = sample_date
  
)

ACI_summary_tree <- ACI_summary[-1:-14,]
#Only want rows with unique ID and Csetpoint
ACI_summary_tree <- distinct(ACI_summary_tree,ID, .keep_all = TRUE)

ACI_sample_date <- unique(subset(ACI_summary_tree, select = c("ID","sample_date")))

#merge meta data
LC_meta <- read.csv(file = "LC_2023/2023_growth_inventory_analysis/LC_9_20_growth_data_cleaned.csv")
ACI_meta_data <- subset(LC_meta, select = c("ID","event_short","construct","construct2","block"))
ACI_meta_data <- inner_join(ACI_meta_data, ACI_sample_date, by = "ID")

ACI_summary_tree <- inner_join(ACI_summary_tree, ACI_meta_data, by=c('ID'))
ACI_summary_tree <- ACI_summary_tree %>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top",
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D"~"control",
  event_short == "CT3"~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))

##create datasheet with all fluoresence dat included
ACI_summary_tree_full_fluor <- na.omit(ACI_summary_tree)
#output
write.csv(ACI_summary_tree,file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_summary_tree_cleaned.csv")
write.csv(ACI_summary_tree_full_fluor, file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_summary_tree_fluor_clean.csv")





#####Sumarize data by event

qt(0.975,8)
ACI_event_summary <- ACI_summary_tree %>% group_by(event_short,Csetpoint) %>% summarise(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  A_qt = qt(0.975,n),
  A_upper = A+(A_qt*A_se),
  A_lower = A-(A_qt*A_se),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = gsw_sd/(sqrt(n)),
  gsw_qt = qt(0.975,n),
  gsw_upper = gsw+(gsw_qt*gsw_se),
  gsw_lower = gsw-(gsw_qt*gsw_se),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = Ci_sd/sqrt(n),
  Ci_qt = qt(0.975,n),
  Ci_upper = Ci+(Ci_qt*Ci_se),
  Ci_lower = Ci-(Ci_qt*Ci_se),
)


ACI_event_summary_full_fluor <- ACI_summary_tree_full_fluor %>% group_by(event_short,Csetpoint) %>% summarise(
  n = n(),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = PhiPS2_sd/sqrt(n),
  PhiPS2_qt = qt(0.975,n),
  PhiPS2_upper = PhiPS2+(PhiPS2_qt*PhiPS2_se),
  PhiPS2_lower = PhiPS2-(PhiPS2_qt*PhiPS2_se),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = ETR_sd/sqrt(n),
  ETR_qt = qt(0.975,n),
  ETR_upper = ETR+(ETR_qt*ETR_se),
  ETR_lower = ETR-(ETR_qt*ETR_se)
)

ACI_event_summary <- inner_join(ACI_event_summary, ACI_event_summary_full_fluor, by = c("event_short","Csetpoint"))
  
ACI_event_summary <- ACI_event_summary%>% mutate(tier = case_when(
  event_short == "5A" | event_short == "5C" | event_short == "4A" ~ "top", 
  event_short == "13-15E" | event_short == "2H" ~ "poor",
  event_short == "16-20" | event_short == "8-9D"~"control",
  event_short == "CT3"~ "WT",
  event_short == "13-15B" | event_short == "1" | event_short == "1C" | event_short == "7"| event_short == "4B"| event_short == "5" ~ "unanalyzed"))

write.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_event_summary.csv",ACI_event_summary)

#### summary of tiers

ACI_tier_summary <- ACI_summary_tree %>% group_by(tier,Csetpoint) %>% summarise(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  A_qt = qt(0.975,n),
  A_upper = A+(A_qt*A_se),
  A_lower = A-(A_qt*A_se),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = gsw_sd/(sqrt(n)),
  gsw_qt = qt(0.975,n),
  gsw_upper = gsw+(gsw_qt*gsw_se),
  gsw_lower = gsw-(gsw_qt*gsw_se),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = Ci_sd/sqrt(n),
  Ci_qt = qt(0.975,n),
  Ci_upper = Ci+(Ci_qt*Ci_se),
  Ci_lower = Ci-(Ci_qt*Ci_se),
)


ACI_tier_summary_full_fluor <- ACI_summary_tree_full_fluor %>% group_by(tier,Csetpoint) %>% summarise(
  n = n(),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = PhiPS2_sd/sqrt(n),
  PhiPS2_qt = qt(0.975,n),
  PhiPS2_upper = PhiPS2+(PhiPS2_qt*PhiPS2_se),
  PhiPS2_lower = PhiPS2-(PhiPS2_qt*PhiPS2_se),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = ETR_sd/sqrt(n),
  ETR_qt = qt(0.975,n),
  ETR_upper = ETR+(ETR_qt*ETR_se),
  ETR_lower = ETR-(ETR_qt*ETR_se)
)

ACI_tier_summary <- inner_join(ACI_tier_summary, ACI_tier_summary_full_fluor, by = c("tier","Csetpoint"))

write.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_tier_summary.csv",ACI_tier_summary)

##### summarize by sample date


ACI_tier_summary_date <- ACI_summary_tree %>% group_by(tier,Csetpoint,sample_date.x) %>% summarise(
  n = n(),
  A_sd = sd(A),
  A = mean(A),
  A_se = A_sd/(sqrt(n)),
  A_qt = qt(0.975,n),
  A_upper = A+(A_qt*A_se),
  A_lower = A-(A_qt*A_se),
  gsw_sd = sd(gsw),
  gsw = mean(gsw),
  gsw_se = gsw_sd/(sqrt(n)),
  gsw_qt = qt(0.975,n),
  gsw_upper = gsw+(gsw_qt*gsw_se),
  gsw_lower = gsw-(gsw_qt*gsw_se),
  Ci_sd = sd(Ci),
  Ci = mean(Ci),
  Ci_se = Ci_sd/sqrt(n),
  Ci_qt = qt(0.975,n),
  Ci_upper = Ci+(Ci_qt*Ci_se),
  Ci_lower = Ci-(Ci_qt*Ci_se),
)


ACI_tier_summary_full_fluor_date <- ACI_summary_tree_full_fluor %>% group_by(tier,Csetpoint,sample_date.x) %>% summarise(
  n = n(),
  PhiPS2_sd = sd(PhiPS2),
  PhiPS2 = mean(PhiPS2),
  PhiPS2_se = PhiPS2_sd/sqrt(n),
  PhiPS2_qt = qt(0.975,n),
  PhiPS2_upper = PhiPS2+(PhiPS2_qt*PhiPS2_se),
  PhiPS2_lower = PhiPS2-(PhiPS2_qt*PhiPS2_se),
  ETR_sd = sd(ETR),
  ETR = mean(ETR),
  ETR_se = ETR_sd/sqrt(n),
  ETR_qt = qt(0.975,n),
  ETR_upper = ETR+(ETR_qt*ETR_se),
  ETR_lower = ETR-(ETR_qt*ETR_se)
)

ACI_tier_summary_date <- inner_join(ACI_tier_summary_date, ACI_tier_summary_full_fluor_date, by = c("tier","Csetpoint"))

write.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/ACI_tier_summary_date.csv",ACI_tier_summary_date)
