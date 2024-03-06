#process fits for A-ci curves and output them in a csv

#install.packages("plantecophys")
#library(plantecophys)
#install.packages("nlstools")
library(nlstools)
library(devtools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(broom)
#install.packages("openxlsx")
library(openxlsx)
library(mgcv)
#install.packages("photosynthesis")
library(photosynthesis)
#read in compiled datasets

AQ_compiled <- read.csv(file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/compiled/light_response_curves_compiled.csv")

#long format

#all_trees_wide <- pivot_wider(ACI_compiled,id_cols = c("ID"), names_from = Csetpoint, values_from = c(A,gsw,Ci))

#fx for viewing data
viewdat = function(dat){plot(y = dat$A, x = dat$Q);abline(h = seq(0,20,1));dat}

#fx for plotting fits

plotfits = function(fit){
  jpeg(filename = paste('',deparse(substitute(fit)),'.jpeg',sep=''),width = 10,height=7,units='in',res = 100)
  plot(fit,main =deparse(substitute(fit)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel),labels = paste('Vcmax =',round(fit$pars[1,1],2)))
  text(pos = 4,x = 0, y = max(fit$df$Amodel)-1,labels = paste('Jmax =',round(fit$pars[2,1],2)))
  dev.off()
}

#options("device")
#create lists to store data in

tree_list = list()
k_sat_list = list()
phi_J_list = list()
theta_J_list = list()
Rd_list = list()
LCP_list = list()
Sample_date = list()

##look at fits, model and graph them and get outputs and assign to storage dataframe

LCOR307 = subset(AQ_compiled, ID == "LCOR-307")
LCOR307_1 = LCOR307[1:8,]

fit_307_1 <- fit_photosynthesis(LCOR307_1, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_307_1 <- coef(fit_307_1) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_307_1) ^ 2)

b_307_1 = coef(fit_307_1)

df_predict_307_1 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_307_1["k_sat"],
      b_307_1["phi_J"],
      b_307_1["theta_J"]
    ) - b_307_1["Rd"]
  )

LCOR307_1_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_307_1) +
  geom_point(data = filter(LCOR307_1)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR307_1_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-307')
b_307_1[1]
k_sat_list <- append(k_sat_list, b_307_1[1])
phi_J_list <- append(phi_J_list, b_307_1[2])
theta_J_list <- append(theta_J_list, b_307_1[3])
Rd_list <- append(Rd_list, b_307_1[4])
LCP_list <- append(LCP_list, LCP_307_1)

###############

LCOR307_2 = LCOR307[9:15,]

fit_307_2 <- fit_photosynthesis(LCOR307_2, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))



#calculate light compensation point
LCP_307_2 <- coef(fit_307_2) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_307_2) ^ 2)

b_307_2 = coef(fit_307_2)

df_predict_307_2 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_307_2["k_sat"],
      b_307_2["phi_J"],
      b_307_2["theta_J"]
    ) - b_307_2["Rd"]
  )

LCOR307_2_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_307_2) +
  geom_point(data = filter(LCOR307_2)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR307_2_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-307_2')
b_307_2[1]
k_sat_list <- append(k_sat_list, b_307_2[1])
phi_J_list <- append(phi_J_list, b_307_2[2])
theta_J_list <- append(theta_J_list, b_307_2[3])
Rd_list <- append(Rd_list, b_307_2[4])
LCP_list <- append(LCP_list, LCP_307_2)

##############

LCOR318 = subset(AQ_compiled, ID == "LCOR-318")

fit_318 <- fit_photosynthesis(LCOR318, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_318 <- coef(fit_318) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_318) ^ 2)

b_318 = coef(fit_318)

df_predict_318 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_318["k_sat"],
      b_318["phi_J"],
      b_318["theta_J"]
    ) - b_318["Rd"]
  )

LCOR318_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_318) +
  geom_point(data = filter(LCOR318)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR318_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-318')
b_318[1]
k_sat_list <- append(k_sat_list, b_318[1])
phi_J_list <- append(phi_J_list, b_318[2])
theta_J_list <- append(theta_J_list, b_318[3])
Rd_list <- append(Rd_list, b_318[4])
LCP_list <- append(LCP_list, LCP_318)

##############

LCOR412 = subset(AQ_compiled, ID == "LCOR-412")

fit_412 <- fit_photosynthesis(LCOR412, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_412 <- coef(fit_412) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_412) ^ 2)

b_412 = coef(fit_412)

df_predict_412 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_412["k_sat"],
      b_412["phi_J"],
      b_412["theta_J"]
    ) - b_412["Rd"]
  )

LCOR412_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_412) +
  geom_point(data = filter(LCOR412)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR412_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-412')
b_412[1]
k_sat_list <- append(k_sat_list, b_412[1])
phi_J_list <- append(phi_J_list, b_412[2])
theta_J_list <- append(theta_J_list, b_412[3])
Rd_list <- append(Rd_list, b_412[4])
LCP_list <- append(LCP_list, LCP_412)

##############

LCOR063 = subset(AQ_compiled, ID == "LCOR-063")

fit_063 <- fit_photosynthesis(LCOR063, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_063 <- coef(fit_063) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_063) ^ 2)

b_063 = coef(fit_063)

df_predict_063 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_063["k_sat"],
      b_063["phi_J"],
      b_063["theta_J"]
    ) - b_063["Rd"]
  )

LCOR063_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_063) +
  geom_point(data = filter(LCOR063)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR063_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-063')
b_063[1]
k_sat_list <- append(k_sat_list, b_063[1])
phi_J_list <- append(phi_J_list, b_063[2])
theta_J_list <- append(theta_J_list, b_063[3])
Rd_list <- append(Rd_list, b_063[4])
LCP_list <- append(LCP_list, LCP_063)

##############

LCOR546 = subset(AQ_compiled, ID == "LCOR-546")

fit_546 <- fit_photosynthesis(LCOR546, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_546 <- coef(fit_546) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_546) ^ 2)

b_546 = coef(fit_546)

df_predict_546 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_546["k_sat"],
      b_546["phi_J"],
      b_546["theta_J"]
    ) - b_546["Rd"]
  )

LCOR546_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_546) +
  geom_point(data = filter(LCOR546)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR546_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-546')
b_546[1]
k_sat_list <- append(k_sat_list, b_546[1])
phi_J_list <- append(phi_J_list, b_546[2])
theta_J_list <- append(theta_J_list, b_546[3])
Rd_list <- append(Rd_list, b_546[4])
LCP_list <- append(LCP_list, LCP_546)

##############

LCOR347 = subset(AQ_compiled, ID == "LCOR-347")

fit_347 <- fit_photosynthesis(LCOR347, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_347 <- coef(fit_347) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_347) ^ 2)

b_347 = coef(fit_347)

df_predict_347 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_347["k_sat"],
      b_347["phi_J"],
      b_347["theta_J"]
    ) - b_347["Rd"]
  )

LCOR347_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_347) +
  geom_point(data = filter(LCOR347)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR347_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-347')
b_347[1]
k_sat_list <- append(k_sat_list, b_347[1])
phi_J_list <- append(phi_J_list, b_347[2])
theta_J_list <- append(theta_J_list, b_347[3])
Rd_list <- append(Rd_list, b_347[4])
LCP_list <- append(LCP_list, LCP_347)

##############

LCOR268 = subset(AQ_compiled, ID == "LCOR-268")

fit_268 <- fit_photosynthesis(LCOR268, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_268 <- coef(fit_268) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_268) ^ 2)

b_268 = coef(fit_268)

df_predict_268 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_268["k_sat"],
      b_268["phi_J"],
      b_268["theta_J"]
    ) - b_268["Rd"]
  )

LCOR268_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_268) +
  geom_point(data = filter(LCOR268)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR268_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-268')
b_268[1]
k_sat_list <- append(k_sat_list, b_268[1])
phi_J_list <- append(phi_J_list, b_268[2])
theta_J_list <- append(theta_J_list, b_268[3])
Rd_list <- append(Rd_list, b_268[4])
LCP_list <- append(LCP_list, LCP_268)

##############

LCOR155 = subset(AQ_compiled, ID == "LCOR-155")

fit_155 <- fit_photosynthesis(LCOR155, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_155 <- coef(fit_155) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_155) ^ 2)

b_155 = coef(fit_155)

df_predict_155 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_155["k_sat"],
      b_155["phi_J"],
      b_155["theta_J"]
    ) - b_155["Rd"]
  )

LCOR155_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_155) +
  geom_point(data = filter(LCOR155)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR155_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-155')
b_155[1]
k_sat_list <- append(k_sat_list, b_155[1])
phi_J_list <- append(phi_J_list, b_155[2])
theta_J_list <- append(theta_J_list, b_155[3])
Rd_list <- append(Rd_list, b_155[4])
LCP_list <- append(LCP_list, LCP_155)

##############

LCOR322 = subset(AQ_compiled, ID == "LCOR-322")

fit_322 <- fit_photosynthesis(LCOR322, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_322 <- coef(fit_322) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_322) ^ 2)

b_322 = coef(fit_322)

df_predict_322 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_322["k_sat"],
      b_322["phi_J"],
      b_322["theta_J"]
    ) - b_322["Rd"]
  )

LCOR322_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_322) +
  geom_point(data = filter(LCOR322)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR322_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-322')
b_322[1]
k_sat_list <- append(k_sat_list, b_322[1])
phi_J_list <- append(phi_J_list, b_322[2])
theta_J_list <- append(theta_J_list, b_322[3])
Rd_list <- append(Rd_list, b_322[4])
LCP_list <- append(LCP_list, LCP_322)

##############

LCOR073 = subset(AQ_compiled, ID == "LCOR-073")

fit_073 <- fit_photosynthesis(LCOR073, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_073 <- coef(fit_073) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_073) ^ 2)

b_073 = coef(fit_073)

df_predict_073 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_073["k_sat"],
      b_073["phi_J"],
      b_073["theta_J"]
    ) - b_073["Rd"]
  )

LCOR073_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_073) +
  geom_point(data = filter(LCOR073)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR073_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-073')
b_073[1]
k_sat_list <- append(k_sat_list, b_073[1])
phi_J_list <- append(phi_J_list, b_073[2])
theta_J_list <- append(theta_J_list, b_073[3])
Rd_list <- append(Rd_list, b_073[4])
LCP_list <- append(LCP_list, LCP_073)

##############

LCOR164 = subset(AQ_compiled, ID == "LCOR-164")

fit_164 <- fit_photosynthesis(LCOR164, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_164 <- coef(fit_164) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_164) ^ 2)

b_164 = coef(fit_164)

df_predict_164 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_164["k_sat"],
      b_164["phi_J"],
      b_164["theta_J"]
    ) - b_164["Rd"]
  )

LCOR164_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_164) +
  geom_point(data = filter(LCOR164)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR164_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-164')
b_164[1]
k_sat_list <- append(k_sat_list, b_164[1])
phi_J_list <- append(phi_J_list, b_164[2])
theta_J_list <- append(theta_J_list, b_164[3])
Rd_list <- append(Rd_list, b_164[4])
LCP_list <- append(LCP_list, LCP_164)

##############

LCOR199 = subset(AQ_compiled, ID == "LCOR-199")

fit_199 <- fit_photosynthesis(LCOR199, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_199 <- coef(fit_199) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_199) ^ 2)

b_199 = coef(fit_199)

df_predict_199 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_199["k_sat"],
      b_199["phi_J"],
      b_199["theta_J"]
    ) - b_199["Rd"]
  )

LCOR199_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_199) +
  geom_point(data = filter(LCOR199)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR199_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-199')
b_199[1]
k_sat_list <- append(k_sat_list, b_199[1])
phi_J_list <- append(phi_J_list, b_199[2])
theta_J_list <- append(theta_J_list, b_199[3])
Rd_list <- append(Rd_list, b_199[4])
LCP_list <- append(LCP_list, LCP_199)

##############

LCOR573 = subset(AQ_compiled, ID == "LCOR-573")

fit_573 <- fit_photosynthesis(LCOR573, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_573 <- coef(fit_573) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_573) ^ 2)

b_573 = coef(fit_573)

df_predict_573 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_573["k_sat"],
      b_573["phi_J"],
      b_573["theta_J"]
    ) - b_573["Rd"]
  )

LCOR573_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_573) +
  geom_point(data = filter(LCOR573)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR573_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-573')
b_573[1]
k_sat_list <- append(k_sat_list, b_573[1])
phi_J_list <- append(phi_J_list, b_573[2])
theta_J_list <- append(theta_J_list, b_573[3])
Rd_list <- append(Rd_list, b_573[4])
LCP_list <- append(LCP_list, LCP_573)

##############

LCOR456 = subset(AQ_compiled, ID == "LCOR-456")

fit_456 <- fit_photosynthesis(LCOR456, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_456 <- coef(fit_456) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_456) ^ 2)

b_456 = coef(fit_456)

df_predict_456 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_456["k_sat"],
      b_456["phi_J"],
      b_456["theta_J"]
    ) - b_456["Rd"]
  )

LCOR456_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_456) +
  geom_point(data = filter(LCOR456)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR456_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-456')
b_456[1]
k_sat_list <- append(k_sat_list, b_456[1])
phi_J_list <- append(phi_J_list, b_456[2])
theta_J_list <- append(theta_J_list, b_456[3])
Rd_list <- append(Rd_list, b_456[4])
LCP_list <- append(LCP_list, LCP_456)

##############

LCOR089 = subset(AQ_compiled, ID == "LCOR-089")

fit_089 <- fit_photosynthesis(LCOR089, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_089 <- coef(fit_089) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_089) ^ 2)

b_089 = coef(fit_089)

df_predict_089 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_089["k_sat"],
      b_089["phi_J"],
      b_089["theta_J"]
    ) - b_089["Rd"]
  )

LCOR089_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_089) +
  geom_point(data = filter(LCOR089)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR089_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-089')
b_089[1]
k_sat_list <- append(k_sat_list, b_089[1])
phi_J_list <- append(phi_J_list, b_089[2])
theta_J_list <- append(theta_J_list, b_089[3])
Rd_list <- append(Rd_list, b_089[4])
LCP_list <- append(LCP_list, LCP_089)

#possibly an outlier, super low. Wasnt included in Aci


##############

LCOR228 = subset(AQ_compiled, ID == "LCOR-228")

fit_228 <- fit_photosynthesis(LCOR228, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_228 <- coef(fit_228) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_228) ^ 2)

b_228 = coef(fit_228)

df_predict_228 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_228["k_sat"],
      b_228["phi_J"],
      b_228["theta_J"]
    ) - b_228["Rd"]
  )

LCOR228_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_228) +
  geom_point(data = filter(LCOR228)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR228_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-228')
b_228[1]
k_sat_list <- append(k_sat_list, b_228[1])
phi_J_list <- append(phi_J_list, b_228[2])
theta_J_list <- append(theta_J_list, b_228[3])
Rd_list <- append(Rd_list, b_228[4])
LCP_list <- append(LCP_list, LCP_228)

##############

LCOR577 = subset(AQ_compiled, ID == "LCOR-577")

fit_577 <- fit_photosynthesis(LCOR577, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_577 <- coef(fit_577) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_577) ^ 2)

b_577 = coef(fit_577)

df_predict_577 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_577["k_sat"],
      b_577["phi_J"],
      b_577["theta_J"]
    ) - b_577["Rd"]
  )

LCOR577_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_577) +
  geom_point(data = filter(LCOR577)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR577_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-577')
b_577[1]
k_sat_list <- append(k_sat_list, b_577[1])
phi_J_list <- append(phi_J_list, b_577[2])
theta_J_list <- append(theta_J_list, b_577[3])
Rd_list <- append(Rd_list, b_577[4])
LCP_list <- append(LCP_list, LCP_577)

##############

LCOR611 = subset(AQ_compiled, ID == "LCOR-611")

fit_611 <- fit_photosynthesis(LCOR611, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_611 <- coef(fit_611) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_611) ^ 2)

b_611 = coef(fit_611)

df_predict_611 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_611["k_sat"],
      b_611["phi_J"],
      b_611["theta_J"]
    ) - b_611["Rd"]
  )

LCOR611_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_611) +
  geom_point(data = filter(LCOR611)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR611_plot

#did not fit well

#save coefficients
#tree_list <- append(tree_list, 'LCOR-611')
#b_611[1]
#k_sat_list <- append(k_sat_list, b_611[1])
#phi_J_list <- append(phi_J_list, b_611[2])
#theta_J_list <- append(theta_J_list, b_611[3])
#Rd_list <- append(Rd_list, b_611[4])
#LCP_list <- append(LCP_list, LCP_611)

##############

LCOR280 = subset(AQ_compiled, ID == "LCOR-280")

fit_280 <- fit_photosynthesis(LCOR280, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_280 <- coef(fit_280) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_280) ^ 2)

b_280 = coef(fit_280)

df_predict_280 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_280["k_sat"],
      b_280["phi_J"],
      b_280["theta_J"]
    ) - b_280["Rd"]
  )

LCOR280_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_280) +
  geom_point(data = filter(LCOR280)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR280_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-280')
b_280[1]
k_sat_list <- append(k_sat_list, b_280[1])
phi_J_list <- append(phi_J_list, b_280[2])
theta_J_list <- append(theta_J_list, b_280[3])
Rd_list <- append(Rd_list, b_280[4])
LCP_list <- append(LCP_list, LCP_280)


##############


LCOR291_2 <- LCOR291[c(8:14),]

fit_291_2 <- fit_photosynthesis(LCOR291_2, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_291_2 <- coef(fit_291_2) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_291_2) ^ 2)

b_291_2 = coef(fit_291_2)

df_predict_291_2 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_291_2["k_sat"],
      b_291_2["phi_J"],
      b_291_2["theta_J"]
    ) - b_291_2["Rd"]
  )

LCOR291_plot_2 <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_291_2) +
  geom_point(data = filter(LCOR291_2)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR291_plot_2

#save coefficients
tree_list <- append(tree_list, 'LCOR-291')
b_291_2[1]
k_sat_list <- append(k_sat_list, b_291_2[1])
phi_J_list <- append(phi_J_list, b_291_2[2])
theta_J_list <- append(theta_J_list, b_291_2[3])
Rd_list <- append(Rd_list, b_291_2[4])
LCP_list <- append(LCP_list, LCP_291_2)

##############

LCOR506 = subset(AQ_compiled, ID == "LCOR-506")

fit_506 <- fit_photosynthesis(LCOR506, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_506 <- coef(fit_506) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_506) ^ 2)

b_506 = coef(fit_506)

df_predict_506 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_506["k_sat"],
      b_506["phi_J"],
      b_506["theta_J"]
    ) - b_506["Rd"]
  )

LCOR506_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_506) +
  geom_point(data = filter(LCOR506)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR506_plot

#poor fits
#save coefficients
#tree_list <- append(tree_list, 'LCOR-506')
#b_506[1]
#k_sat_list <- append(k_sat_list, b_506[1])
#phi_J_list <- append(phi_J_list, b_506[2])
#theta_J_list <- append(theta_J_list, b_506[3])
#Rd_list <- append(Rd_list, b_506[4])
#LCP_list <- append(LCP_list, LCP_506)

##############

LCOR501 = subset(AQ_compiled, ID == "LCOR-501")

fit_501 <- fit_photosynthesis(LCOR501, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_501 <- coef(fit_501) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_501) ^ 2)

b_501 = coef(fit_501)

df_predict_501 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_501["k_sat"],
      b_501["phi_J"],
      b_501["theta_J"]
    ) - b_501["Rd"]
  )

LCOR501_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_501) +
  geom_point(data = filter(LCOR501)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR501_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-501')
b_501[1]
k_sat_list <- append(k_sat_list, b_501[1])
phi_J_list <- append(phi_J_list, b_501[2])
theta_J_list <- append(theta_J_list, b_501[3])
Rd_list <- append(Rd_list, b_501[4])
LCP_list <- append(LCP_list, LCP_501)

##############

LCOR211 = subset(AQ_compiled, ID == "LCOR-211")

fit_211 <- fit_photosynthesis(LCOR211, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_211 <- coef(fit_211) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_211) ^ 2)

b_211 = coef(fit_211)

df_predict_211 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_211["k_sat"],
      b_211["phi_J"],
      b_211["theta_J"]
    ) - b_211["Rd"]
  )

LCOR211_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_211) +
  geom_point(data = filter(LCOR211)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR211_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-211')
b_211[1]
k_sat_list <- append(k_sat_list, b_211[1])
phi_J_list <- append(phi_J_list, b_211[2])
theta_J_list <- append(theta_J_list, b_211[3])
Rd_list <- append(Rd_list, b_211[4])
LCP_list <- append(LCP_list, LCP_211)

##############

LCOR274 = subset(AQ_compiled, ID == "LCOR-274")

fit_274 <- fit_photosynthesis(LCOR274, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_274 <- coef(fit_274) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_274) ^ 2)

b_274 = coef(fit_274)

df_predict_274 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_274["k_sat"],
      b_274["phi_J"],
      b_274["theta_J"]
    ) - b_274["Rd"]
  )

LCOR274_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_274) +
  geom_point(data = filter(LCOR274)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR274_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-274')
b_274[1]
k_sat_list <- append(k_sat_list, b_274[1])
phi_J_list <- append(phi_J_list, b_274[2])
theta_J_list <- append(theta_J_list, b_274[3])
Rd_list <- append(Rd_list, b_274[4])
LCP_list <- append(LCP_list, LCP_274)

##############

LCOR231 = subset(AQ_compiled, ID == "LCOR-231")

fit_231 <- fit_photosynthesis(LCOR231, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_231 <- coef(fit_231) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_231) ^ 2)

b_231 = coef(fit_231)

df_predict_231 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_231["k_sat"],
      b_231["phi_J"],
      b_231["theta_J"]
    ) - b_231["Rd"]
  )

LCOR231_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_231) +
  geom_point(data = filter(LCOR231)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR231_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-231')
b_231[1]
k_sat_list <- append(k_sat_list, b_231[1])
phi_J_list <- append(phi_J_list, b_231[2])
theta_J_list <- append(theta_J_list, b_231[3])
Rd_list <- append(Rd_list, b_231[4])
LCP_list <- append(LCP_list, LCP_231)

##############

LCOR101 = subset(AQ_compiled, ID == "LCOR-101")

fit_101 <- fit_photosynthesis(LCOR101, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_101 <- coef(fit_101) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_101) ^ 2)

b_101 = coef(fit_101)

df_predict_101 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_101["k_sat"],
      b_101["phi_J"],
      b_101["theta_J"]
    ) - b_101["Rd"]
  )

LCOR101_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_101) +
  geom_point(data = filter(LCOR101)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR101_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-101')
b_101[1]
k_sat_list <- append(k_sat_list, b_101[1])
phi_J_list <- append(phi_J_list, b_101[2])
theta_J_list <- append(theta_J_list, b_101[3])
Rd_list <- append(Rd_list, b_101[4])
LCP_list <- append(LCP_list, LCP_101)

##############

LCOR204 = subset(AQ_compiled, ID == "LCOR-204")

fit_204 <- fit_photosynthesis(LCOR204, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_204 <- coef(fit_204) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_204) ^ 2)

b_204 = coef(fit_204)

df_predict_204 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_204["k_sat"],
      b_204["phi_J"],
      b_204["theta_J"]
    ) - b_204["Rd"]
  )

LCOR204_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_204) +
  geom_point(data = filter(LCOR204)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR204_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-204')
b_204[1]
k_sat_list <- append(k_sat_list, b_204[1])
phi_J_list <- append(phi_J_list, b_204[2])
theta_J_list <- append(theta_J_list, b_204[3])
Rd_list <- append(Rd_list, b_204[4])
LCP_list <- append(LCP_list, LCP_204)

##############

LCOR203 = subset(AQ_compiled, ID == "LCOR-203")

fit_203 <- fit_photosynthesis(LCOR203, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_203 <- coef(fit_203) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_203) ^ 2)

b_203 = coef(fit_203)

df_predict_203 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_203["k_sat"],
      b_203["phi_J"],
      b_203["theta_J"]
    ) - b_203["Rd"]
  )

LCOR203_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_203) +
  geom_point(data = filter(LCOR203)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR203_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-203')
b_203[1]
k_sat_list <- append(k_sat_list, b_203[1])
phi_J_list <- append(phi_J_list, b_203[2])
theta_J_list <- append(theta_J_list, b_203[3])
Rd_list <- append(Rd_list, b_203[4])
LCP_list <- append(LCP_list, LCP_203)

##############

LCOR610 = subset(AQ_compiled, ID == "LCOR-610")

fit_610 <- fit_photosynthesis(LCOR610, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_610 <- coef(fit_610) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_610) ^ 2)

b_610 = coef(fit_610)

df_predict_610 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_610["k_sat"],
      b_610["phi_J"],
      b_610["theta_J"]
    ) - b_610["Rd"]
  )

LCOR610_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_610) +
  geom_point(data = filter(LCOR610)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR610_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-610')
b_610[1]
k_sat_list <- append(k_sat_list, b_610[1])
phi_J_list <- append(phi_J_list, b_610[2])
theta_J_list <- append(theta_J_list, b_610[3])
Rd_list <- append(Rd_list, b_610[4])
LCP_list <- append(LCP_list, LCP_610)

##############

LCOR275 = subset(AQ_compiled, ID == "LCOR-275")

fit_275 <- fit_photosynthesis(LCOR275, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_275 <- coef(fit_275) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_275) ^ 2)

b_275 = coef(fit_275)

df_predict_275 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_275["k_sat"],
      b_275["phi_J"],
      b_275["theta_J"]
    ) - b_275["Rd"]
  )

LCOR275_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_275) +
  geom_point(data = filter(LCOR275)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR275_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-275')
b_275[1]
k_sat_list <- append(k_sat_list, b_275[1])
phi_J_list <- append(phi_J_list, b_275[2])
theta_J_list <- append(theta_J_list, b_275[3])
Rd_list <- append(Rd_list, b_275[4])
LCP_list <- append(LCP_list, LCP_275)

##############

LCOR516 = subset(AQ_compiled, ID == "LCOR-516")

fit_516 <- fit_photosynthesis(LCOR516, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_516 <- coef(fit_516) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_516) ^ 2)

b_516 = coef(fit_516)

df_predict_516 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_516["k_sat"],
      b_516["phi_J"],
      b_516["theta_J"]
    ) - b_516["Rd"]
  )

LCOR516_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_516) +
  geom_point(data = filter(LCOR516)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR516_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-516')
b_516[1]
k_sat_list <- append(k_sat_list, b_516[1])
phi_J_list <- append(phi_J_list, b_516[2])
theta_J_list <- append(theta_J_list, b_516[3])
Rd_list <- append(Rd_list, b_516[4])
LCP_list <- append(LCP_list, LCP_516)

##############

LCOR003 = subset(AQ_compiled, ID == "LCOR-003")

fit_003 <- fit_photosynthesis(LCOR003, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_003 <- coef(fit_003) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_003) ^ 2)

b_003 = coef(fit_003)

df_predict_003 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_003["k_sat"],
      b_003["phi_J"],
      b_003["theta_J"]
    ) - b_003["Rd"]
  )

LCOR003_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_003) +
  geom_point(data = filter(LCOR003)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR003_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-003')
b_003[1]
k_sat_list <- append(k_sat_list, b_003[1])
phi_J_list <- append(phi_J_list, b_003[2])
theta_J_list <- append(theta_J_list, b_003[3])
Rd_list <- append(Rd_list, b_003[4])
LCP_list <- append(LCP_list, LCP_003)

##############

LCOR578 = subset(AQ_compiled, ID == "LCOR-578")

fit_578 <- fit_photosynthesis(LCOR578, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_578 <- coef(fit_578) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_578) ^ 2)

b_578 = coef(fit_578)

df_predict_578 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_578["k_sat"],
      b_578["phi_J"],
      b_578["theta_J"]
    ) - b_578["Rd"]
  )

LCOR578_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_578) +
  geom_point(data = filter(LCOR578)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR578_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-578')
b_578[1]
k_sat_list <- append(k_sat_list, b_578[1])
phi_J_list <- append(phi_J_list, b_578[2])
theta_J_list <- append(theta_J_list, b_578[3])
Rd_list <- append(Rd_list, b_578[4])
LCP_list <- append(LCP_list, LCP_578)

##############

LCOR417 = subset(AQ_compiled, ID == "LCOR-417")

fit_417 <- fit_photosynthesis(LCOR417, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_417 <- coef(fit_417) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_417) ^ 2)

b_417 = coef(fit_417)

df_predict_417 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_417["k_sat"],
      b_417["phi_J"],
      b_417["theta_J"]
    ) - b_417["Rd"]
  )

LCOR417_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_417) +
  geom_point(data = filter(LCOR417)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR417_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-417')
b_417[1]
k_sat_list <- append(k_sat_list, b_417[1])
phi_J_list <- append(phi_J_list, b_417[2])
theta_J_list <- append(theta_J_list, b_417[3])
Rd_list <- append(Rd_list, b_417[4])
LCP_list <- append(LCP_list, LCP_417)

##############

LCOR303 = subset(AQ_compiled, ID == "LCOR-303")

fit_303 <- fit_photosynthesis(LCOR303, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_303 <- coef(fit_303) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_303) ^ 2)

b_303 = coef(fit_303)

df_predict_303 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_303["k_sat"],
      b_303["phi_J"],
      b_303["theta_J"]
    ) - b_303["Rd"]
  )

LCOR303_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_303) +
  geom_point(data = filter(LCOR303)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR303_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-303')
b_303[1]
k_sat_list <- append(k_sat_list, b_303[1])
phi_J_list <- append(phi_J_list, b_303[2])
theta_J_list <- append(theta_J_list, b_303[3])
Rd_list <- append(Rd_list, b_303[4])
LCP_list <- append(LCP_list, LCP_303)

##############

LCOR502 = subset(AQ_compiled, ID == "LCOR-502")

#fit_502 <- fit_photosynthesis(LCOR502, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))
#too low to fit

#calculate light compensation point
#LCP_502 <- coef(fit_502) |>
 # t() |>
  #as.data.frame() |>
  #mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  #sum(resid(fit_502) ^ 2)

#b_502 = coef(fit_502)

#df_predict_502 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
 # mutate(
  #  A = marshall_biscoe_1980(
   #   Q_abs = Qin,
    #  k_sat = b_502["k_sat"],
     # b_502["phi_J"],
      #b_502["theta_J"]
    #) - b_502["Rd"]
  #)

#LCOR502_plot <- ggplot(mapping = aes(Qin, A)) +
 # geom_line(data = df_predict_502) +
  #geom_point(data = filter(LCOR502)) +
  #labs(
   # x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    #y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  #) +
  #theme_bw()
#LCOR502_plot

#save coefficients
#tree_list <- append(tree_list, 'LCOR-502')
#b_502[1]
#k_sat_list <- append(k_sat_list, b_502[1])
#phi_J_list <- append(phi_J_list, b_502[2])
#theta_J_list <- append(theta_J_list, b_502[3])
#Rd_list <- append(Rd_list, b_502[4])
#LCP_list <- append(LCP_list, LCP_502)

##############

LCOR314 = subset(AQ_compiled, ID == "LCOR-314")

fit_314 <- fit_photosynthesis(LCOR314, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_314 <- coef(fit_314) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_314) ^ 2)

b_314 = coef(fit_314)

df_predict_314 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_314["k_sat"],
      b_314["phi_J"],
      b_314["theta_J"]
    ) - b_314["Rd"]
  )

LCOR314_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_314) +
  geom_point(data = filter(LCOR314)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR314_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-314')
b_314[1]
k_sat_list <- append(k_sat_list, b_314[1])
phi_J_list <- append(phi_J_list, b_314[2])
theta_J_list <- append(theta_J_list, b_314[3])
Rd_list <- append(Rd_list, b_314[4])
LCP_list <- append(LCP_list, LCP_314)

##############

LCOR092 = subset(AQ_compiled, ID == "LCOR-092")

fit_092 <- fit_photosynthesis(LCOR092, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_092 <- coef(fit_092) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_092) ^ 2)

b_092 = coef(fit_092)

df_predict_092 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_092["k_sat"],
      b_092["phi_J"],
      b_092["theta_J"]
    ) - b_092["Rd"]
  )

LCOR092_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_092) +
  geom_point(data = filter(LCOR092)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR092_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-092')
b_092[1]
k_sat_list <- append(k_sat_list, b_092[1])
phi_J_list <- append(phi_J_list, b_092[2])
theta_J_list <- append(theta_J_list, b_092[3])
Rd_list <- append(Rd_list, b_092[4])
LCP_list <- append(LCP_list, LCP_092)

##############

LCOR065 = subset(AQ_compiled, ID == "LCOR-065")

fit_065 <- fit_photosynthesis(LCOR065, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_065 <- coef(fit_065) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_065) ^ 2)

b_065 = coef(fit_065)

df_predict_065 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_065["k_sat"],
      b_065["phi_J"],
      b_065["theta_J"]
    ) - b_065["Rd"]
  )

LCOR065_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_065) +
  geom_point(data = filter(LCOR065)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR065_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-065')
b_065[1]
k_sat_list <- append(k_sat_list, b_065[1])
phi_J_list <- append(phi_J_list, b_065[2])
theta_J_list <- append(theta_J_list, b_065[3])
Rd_list <- append(Rd_list, b_065[4])
LCP_list <- append(LCP_list, LCP_065)

##############

LCOR165 = subset(AQ_compiled, ID == "LCOR-165")

fit_165 <- fit_photosynthesis(LCOR165, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_165 <- coef(fit_165) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_165) ^ 2)

b_165 = coef(fit_165)

df_predict_165 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_165["k_sat"],
      b_165["phi_J"],
      b_165["theta_J"]
    ) - b_165["Rd"]
  )

LCOR165_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_165) +
  geom_point(data = filter(LCOR165)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR165_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-165')
b_165[1]
k_sat_list <- append(k_sat_list, b_165[1])
phi_J_list <- append(phi_J_list, b_165[2])
theta_J_list <- append(theta_J_list, b_165[3])
Rd_list <- append(Rd_list, b_165[4])
LCP_list <- append(LCP_list, LCP_165)

##############

##############

LCOR040 = subset(AQ_compiled, ID == "LCOR-040")

fit_040 <- fit_photosynthesis(LCOR040, .photo_fun = "aq_response",.vars = list(.A = "A",.Q = "Qin"))


#calculate light compensation point
LCP_040 <- coef(fit_040) |>
  t() |>
  as.data.frame() |>
  mutate(LCP = ((Rd) * (Rd * theta_J - k_sat) / (phi_J * (Rd - k_sat)))) |>
  
  ## Calculate residual sum-of-squares
  sum(resid(fit_040) ^ 2)

b_040 = coef(fit_040)

df_predict_040 = data.frame(Qin = seq(0, 0.84 * 2000, length.out = 100)) |>
  mutate(
    A = marshall_biscoe_1980(
      Q_abs = Qin,
      k_sat = b_040["k_sat"],
      b_040["phi_J"],
      b_040["theta_J"]
    ) - b_040["Rd"]
  )

LCOR040_plot <- ggplot(mapping = aes(Qin, A)) +
  geom_line(data = df_predict_040) +
  geom_point(data = filter(LCOR040)) +
  labs(
    x = expression("Irradiance (" * mu * mol ~ m^{-2} ~ s^{-1} * ")"),
    y = expression(A[net] ~ "(" * mu * mol ~ m^{-2} ~ s^{-1} * ")")
  ) +
  theme_bw()
LCOR040_plot

#save coefficients
tree_list <- append(tree_list, 'LCOR-040')
b_040[1]
k_sat_list <- append(k_sat_list, b_040[1])
phi_J_list <- append(phi_J_list, b_040[2])
theta_J_list <- append(theta_J_list, b_040[3])
Rd_list <- append(Rd_list, b_040[4])
LCP_list <- append(LCP_list, LCP_040)


##########

tree_list <- unlist(tree_list)
tree_list <- as.data.frame(as.matrix(tree_list))
k_sat_list <- unlist(k_sat_list)
k_sat_list <- as.data.frame(as.matrix(k_sat_list))
k_sat_list2 <- filter(k_sat_list, rownames(k_sat_list) %in% rownames(tree_list))
phi_J_list <- unlist(theta_J_list)
phi_J_list <- as.data.frame(as.matrix(phi_J_list))
phi_J_list2 <- filter(phi_J_list, rownames(phi_J_list) %in% rownames(tree_list))
theta_J_list <- unlist(theta_J_list)
theta_J_list <- as.data.frame(as.matrix(theta_J_list))
theta_J_list2 <- filter(theta_J_list, rownames(theta_J_list) %in% rownames(tree_list))
Rd_list <- unlist(Rd_list)
Rd_list <- as.data.frame(as.matrix(Rd_list))
Rd_list2 <- filter(Rd_list, rownames(Rd_list) %in% rownames(tree_list))


AQ_pars <- cbind(tree_list,k_sat_list,phi_J_list,theta_J_list,Rd_list)

colnames(AQ_pars) = c('ID', 'k_sat','phi_J',"theta_J",'Rd')
AQ_pars$ID_2 <- substr(AQ_pars$ID,1,8)
str(efficiency_pars)

AQ_pars2 <- AQ_pars %>% group_by(ID_2) %>% summarize(
  k_sat = mean(k_sat),
  phi_J = mean(phi_J),
  theta_J = mean(theta_J),
  Rd = mean(Rd)
)

AQ_pars2$ID <- AQ_pars2$ID_2



LC_meta <- read.csv("LC_2022/2022_growth&inventory_analylsis/growth_analysis/3_22_growth_cleaned_II.csv")
LC_meta <- subset(LC_meta, select = c("row", "column","ID","event","event_short","block","construct","construct2"))


AQ_parameters_list <- inner_join(AQ_pars2, LC_meta, by = "ID")
sample_date <- subset(AQ_compiled, select = c("ID","sample_date"))

AQ_parameters_list2 <- distinct(inner_join(AQ_parameters_list, sample_date, by = "ID"))


write.csv(AQ_parameters_list2, file = "LC_2023/2023_physiology_analysis/LC_2023_Response_curves/AQ_parameters_list.csv")
