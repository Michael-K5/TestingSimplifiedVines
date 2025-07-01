# Script for performing a quantile regression to test for the simplifying assumption.
source("R/simulate_non_simplified_vine.R")
source("R/classifier_methods.R")
library(ggplot2)
# Parameters to determine, which model and data to load.
last_data_simulation_date <- "2025-06-25"
last_train_date <- "2025-06-25"
data_dim <- 5
nu <- 1 # T_n/T_c, the ratio of noise to observed samples during training
# load the Neural Network Classifier
model_path <- paste0("models/NN_", data_dim, "d_", last_train_date, ".keras")
model <- load_model_hdf5(model_path)
# load the fitted simplified vine
cop_path <- paste0("models/copula_",data_dim,"d_", last_train_date,".rds")
fitted_cop <- readRDS(file = cop_path)
# load the original data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"d_",last_data_simulation_date,".csv")
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames


cor_facs <- correction_factors(model, obs=orig_data, nu=nu)
# how many samples to remove
remove_top <- floor(0.02 * length(cor_facs))
# remove the highest observations, by first sorting
# and then removing the last remove_top observations
cor_facs_no_outliers <- (sort(cor_facs)[1:(length(cor_facs) - remove_top)])
plot(cor_facs_no_outliers, xlab="", ylab="Correction Factors")
cor_facs_hist_KDE_plot <- ggplot(data.frame(x = cor_facs_no_outliers), aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "darkblue", linewidth = 1) +
  labs(title = "Histogram and KDE of the correction factors",
       x = "Correction Factors", y = "Density") +
  theme_minimal()
cor_facs_hist_KDE_plot
ggsave(
  filename = "CorrectionFactorsHistKDE20250625.png",
  plot = cor_facs_hist_KDE_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg="white"
)
log_cor_facs_hist_KDE_plot <- ggplot(data.frame(x = log(cor_facs)), aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "darkblue", linewidth = 1) +
  labs(title = "Histogram and KDE of the correction factors",
       x = "Log Correction Factors", y = "Density") +
  theme_minimal()
log_cor_facs_hist_KDE_plot
ggsave(
  filename = "LogCorrectionFactorsHistKDE20250625.png",
  plot = log_cor_facs_hist_KDE_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg="white"
)
# Noise Densities
noise_densities <- dvinecop(orig_data, fitted_cop)
plot(noise_densities, cor_facs, main="Noise Density vs Correction Factors",
     xlab="Simplified Vine Density (Noise)", ylab="Correction_factors")
log_noise_densities <- log(noise_densities)
log_cor_facs <- log(cor_facs)
plot(log_noise_densities, log_cor_facs, main="Log Noise Density vs Log Correction Factors",
     xlab="Log Simplified Vine Density (Noise)", ylab="Log Correction Factors")
df_temp <- data.frame(log_noise_dens = log_noise_densities, log_cor_factors = log_cor_facs)

hex_plot_log_noise_log_cor_facs <- ggplot(df_temp, aes(x = log_noise_dens, y = log_cor_factors)) +
  geom_hex(binwidth = c(0.5, 0.5)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title= "Log Noise Density vs Log Correction Factors",
    x = "Log Simplified Vine Density (Noise)",
    y = "Log Correction Factors",
    fill = "Count"
  ) +
  theme_minimal()
hex_plot_log_noise_log_cor_facs
ggsave(
  filename = "LogNoiseDensVsLogCorFacs20250625.png",
  plot = hex_plot_log_noise_log_cor_facs,
  width = 8,
  height = 8,
  dpi = 300,
  bg="white"
)

model_densities <- non_param_cop(model=model, fitted_vine=fitted_cop, obs=orig_data, nu=nu)
log_model_densities <- log(model_densities)
df_temp <- data.frame(log_noise_dens = log_noise_densities, log_model_dens = log_model_densities)
noise_vs_model_dens_log_plot <- ggplot(df_temp, aes(x = log_noise_dens, y = log_model_dens)) +
  geom_hex(binwidth = c(0.5, 0.5)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(
    title= "Log Noise Density vs Log Model Density",
    x = "Log Simplified Vine Density (Noise)",
    y = "Log Model Density (NCE)",
    fill = "Count"
  ) +
  theme_minimal()
noise_vs_model_dens_log_plot
ggsave(
  filename = "LogNoiseVsLogModelDens20250625.png",
  plot = noise_vs_model_dens_log_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg="white"
)
# # Try on simplified data
# simp_data <- rvinecop(10000, fitted_cop)
# model_densities_simp <- non_param_cop(model=model, fitted_vine=fitted_cop, obs=simp_data, nu=nu)
# log_model_densities_simp <- log(model_densities_simp)
# log_noise_dens_simp <- log(dvinecop(simp_data, fitted_cop))
# df_temp <- data.frame(log_noise_dens_simp = log_noise_dens_simp, log_model_dens_simp = log_model_densities_simp)
# noise_vs_model_dens_log_plot_simp <- ggplot(df_temp, aes(x = log_noise_dens_simp, y = log_model_dens_simp)) +
#   geom_hex(binwidth = c(0.5, 0.5)) +
#   scale_fill_gradient(low = "lightblue", high = "darkblue") +
#   geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
#   labs(
#     title= "Log Noise vs Log Model Density on Simp Data",
#     x = "Log Simplified Vine Density (Noise)",
#     y = "Log Model Density (NCE)",
#     fill = "Count"
#   ) +
#   theme_minimal()
# noise_vs_model_dens_log_plot_simp

# Norm U against log model and noise densities:
orig_data_norms <- sqrt(rowSums(orig_data^2))
# theo_max <- sqrt(ncol(orig_data))
# theo_max
# max(orig_data_norms)
df_temp <- data.frame(u_norms = orig_data_norms, log_noise_dens = log_noise_densities, log_model_dens = log_model_densities)
p1 <- ggplot(df_temp, aes(x = u_norms, y = log_noise_dens)) +
  geom_hex(binwidth = c(0.05, 0.5)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title= "Log Noise Density vs Log Model Density",
    x = "2-Norm of Data U",
    y = "Log Simplified Vine Density (Noise)",
    fill = "Count"
  ) +
  theme_minimal()
p2 <- ggplot(df_temp, aes(x = u_norms, y = log_model_dens)) +
  geom_hex(binwidth = c(0.05, 0.5)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title= "Log Noise Density vs Log Model Density",
    x = "2-Norm of Data U",
    y = "Log Model Density (NCE)",
    fill = "Count"
  ) +
  theme_minimal()
library(patchwork)
combined_plot <- p1 / p2 # stack plots vertically
combined_plot
ggsave(
  filename = "NormOfDataVsLogNoiseAndModelDens20250625.png",
  plot = combined_plot,
  width = 6,
  height = 10,
  dpi = 300,
  bg="white"
)
# Huks Paper
# Using the Paper of Huk: c/(1+c) = P[u_i is from c versus independent]
simp_vs_ind <- noise_densities / (1+noise_densities)
model_vs_ind <- model_densities / (1+model_densities)
df_temp <- data.frame(simplified_vs_indep = simp_vs_ind, model_vs_indep = model_vs_ind)
simp_ind_vs_model_ind <- ggplot(df_temp, aes(x = simplified_vs_indep, y = model_vs_indep)) +
  geom_hex(binwidth = c(0.03, 0.03)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title= "Plot of c/(1+c) for simplified and NCE model",
    x = "Simplified versus independent",
    y = "Model versus independent",
    fill = "Count"
  ) +
  theme_minimal()
simp_ind_vs_model_ind
ggsave(
  filename = "HukSimpIndVsModelInd20250625.png",
  plot = simp_ind_vs_model_ind,
  width = 10,
  height = 10,
  dpi = 300,
  bg="white"
)
# temp <- non_param_cop(orig_data)
int_val <- compute_integral(model=model,
                            fitted_vine = fitted_cop,
                            nu=nu,
                            data_dim_if_unif = as.integer(data_dim),
                            n_samples=20000,
                            user_info=TRUE)
int_val
# temp_norm <- temp / int_val

# Wilcox and t-test for log_cor_facs = log(c_model) - log(c_noise)
cor_facs <- correction_factors(model, obs=orig_data, nu=nu)
log_cor_facs <- log(cor_facs)
wilcox.test(log_cor_facs, alternative="greater",exact=FALSE,conf.int=TRUE,conf.level=0.95)
t.test(log_cor_facs, alternative="greater", mu=0)
# Equivalent would be the following:
# noise_densities <- dvinecop(orig_data, fitted_cop)
# log_noise_densities <- log(noise_densities)
# model_densities <- non_param_cop(model=model, fitted_vine=fitted_cop, obs=orig_data, nu=nu)
# log_model_densities <- log(model_densities)
# wilcox.test(x=log_model_densities, y=log_noise_densities, alternative="greater",exact=FALSE, paired=TRUE,
#             conf.int=TRUE, conf.level=0.95, mu=0)
# t.test(x=log_model_densities, y=log_noise_densities, alternative="greater", paired=TRUE, cond.level=0.95, mu=0)
# End Wilcox and t-test
# Likelihood ratio test (Vuong)
# evaluate log likelihoods of original data (on train and test set fractions of orig_data)
log_lik_simp_orig <- sum(log(dvinecop(orig_data, fitted_cop)))
log_lik_NN_orig <- sum(log(non_param_cop(
  model=model,
  fitted_vine=fitted_cop,
  obs=orig_data,
  nu=nu)))
# get number of model parameters and degrees of freedom for the log likelihood ratio test
NN_num_params <- count_NN_params(weights=model$weights)
simp_cop_num_params <- get_num_cop_params(fitted_cop)
# NCE model has params from simp_cop and NN, simplified model only from simp_cop
# so the degrees of freedom are just the number of parameters in the neural network
deg_free <- NN_num_params
# log likelihood ratio test statistic
LRT_stat <- -2*(as.numeric(log_lik_simp_orig) - as.numeric(log_lik_NN_orig))
# compute p(chisq(df=deg_free) > LRT_stat) (for train and test set)
p_value_train <- pchisq(LRT_stat, df = deg_free, lower.tail = FALSE)
p_value_train
# compute AIC and BIC for the NN and simplified model
AIC_NN <- 2*(NN_num_params + simp_cop_num_params) - 2* log_lik_NN_orig
AIC_simp <- 2*simp_cop_num_params - 2*log_lik_simp_orig
BIC_NN <- log(nrow(orig_data))*(NN_num_params + simp_cop_num_params) - 2* log_lik_NN_orig
BIC_simp <- log(nrow(orig_data))*simp_cop_num_params - 2*log_lik_simp_orig
AIC_NN
AIC_simp
BIC_NN
BIC_simp
# End likelihood ratio test


# Quantile regression:
# Compute the values G(u, eta) from the thesis
g_vals <- compute_gvals(model,orig_data,nu=nu)
# D-Vine quantile regression:
print("D-Vine quantile regression, nonparametric")
bottom_quantiles <- seq(0.01,0.2,0.01)
top_quantiles <- seq(0.8,0.99,0.01)
output_qreg <- perform_quant_reg(
  orig_data=orig_data,
  g_vals=g_vals,
  family_set_name = "nonparametric",
  bottom_quantile_levels = bottom_quantiles,
  top_quantile_levels = top_quantiles
)

alternative_better_dvine_nonpar <- output_qreg[[1]]
simp_better_dvine_nonpar <- output_qreg[[2]]
dvine_nonpar_loss <- output_qreg[[3]]
# 1 percent, 5 percent, 10 percent etc.
dvine_nonpar_loss[c(1,5,10,15,20,21,26,31,36,40)]
for(i in 1:length(bottom_quantiles)){
  print(paste0("At the alpha = ",
               bottom_quantiles[i],
               " level, the alternative model is better in ",
               alternative_better_dvine_nonpar[i],
               " cases, which is a fraction of ",
               alternative_better_dvine_nonpar[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_quantiles)){
  print(paste0("At the alpha = ",
               1 - top_quantiles[i],
               " level, the simplified model is better in ",
               simp_better_dvine_nonpar[i],
               " cases, which is a fraction of ",
               simp_better_dvine_nonpar[i] / nrow(orig_data),
               "."))
}

# D-Vine quantile regression: onepar
print("D-Vine quantile regression, onepar")
output_qreg <- perform_quant_reg(
  orig_data,
  g_vals,
  family_set_name = "onepar",
  bottom_quantile_levels = bottom_quantiles,
  top_quantile_levels = top_quantiles
)

alternative_better_dvine_onepar <- output_qreg[[1]]
simp_better_dvine_onepar <- output_qreg[[2]]
dvine_onepar_loss <- output_qreg[[3]]
for(i in 1:length(bottom_quantiles)){
  print(paste0("At the alpha = ",
               bottom_quantiles[i],
               " level, the alternative model is better in ",
               alternative_better_dvine_onepar[i],
               " cases, which is a fraction of ",
               alternative_better_dvine_onepar[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_quantiles)){
  print(paste0("At the alpha = ",
               1 - top_quantiles[i],
               " level, the simplified model is better in ",
               simp_better_dvine_onepar[i],
               " cases, which is a fraction of ",
               simp_better_dvine_onepar[i] / nrow(orig_data),
               "."))
}

# D-Vine quantile regression: parametric
print("D-Vine quantile regression, parametric")
output_qreg <- perform_quant_reg(
  orig_data,
  g_vals,
  family_set_name = "parametric",
  bottom_quantile_levels = bottom_quantiles,
  top_quantile_levels = top_quantiles
)

alternative_better_dvine_par <- output_qreg[[1]]
simp_better_dvine_par <- output_qreg[[2]]
dvine_par_loss <- output_qreg[[3]]
for(i in 1:length(bottom_quantiles)){
  print(paste0("At the alpha = ",
               bottom_quantiles[i],
               " level, the alternative model is better in ",
               alternative_better_dvine_par[i],
               " cases, which is a fraction of ",
               alternative_better_dvine_par[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_quantiles)){
  print(paste0("At the alpha = ",
               1 - top_quantiles[i],
               " level, the simplified model is better in ",
               simp_better_dvine_par[i],
               " cases, which is a fraction of ",
               simp_better_dvine_par[i] / nrow(orig_data),
               "."))
}


print("Linear Quantile Regression")
bottom_quantiles_lin <-  seq(0.01,0.2,0.01)
top_quantiles_lin <- seq(0.8,0.99,0.01)
lin_qreg_output <- perform_linear_quant_reg(
  orig_data,
  g_vals,
  bottom_quantiles_lin=bottom_quantiles_lin,
  top_quantiles_lin=top_quantiles_lin)
alternative_better_lin <- lin_qreg_output[[1]]
simp_better_lin <- lin_qreg_output[[2]]
lin_train_loss <- lin_qreg_output[[3]]
lin_test_loss <- lin_qreg_output[[4]]
# loss for 1 percent, 5 percent, 10 percent etc.
lin_train_loss[c(1,5,10,15,20,21,26,31,36,40)]
lin_test_loss[c(1,5,10,15,20,21,26,31,36,40)]
for(i in 1:length(bottom_quantiles_lin)){
  print(paste0("At the alpha = ",
               bottom_quantiles_lin[i],
               " level, the alternative model is better in ",
               alternative_better_lin[i],
               " cases, which is a fraction of ",
               alternative_better_lin[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_quantiles_lin)){
  print(paste0("At the alpha = ",
               1 - top_quantiles_lin[i],
               " level, the simplified model is better in ",
               simp_better_lin[i],
               " cases, which is a fraction of ",
               simp_better_lin[i] / nrow(orig_data),
               "."))
}
if(length(bottom_quantiles_lin) == length(top_quantiles_lin)){
  for(i in 1:length(bottom_quantiles_lin)){
    if(abs(bottom_quantiles_lin[i]- (1 - top_quantiles_lin[length(top_quantiles_lin) - i+ 1])) < 0.0001){
      print(paste0("At the alpha = ",
                   bottom_quantiles_lin[i],
                   " level, the test is inconclusive in ",
                   nrow(orig_data) - (alternative_better_lin[i] + simp_better_lin[length(top_quantiles_lin) - i+ 1]),
                   " cases, which is a fraction of ",
                   (nrow(orig_data) - (alternative_better_lin[i] + simp_better_lin[length(top_quantiles_lin) - i+ 1]))/nrow(orig_data),
                   "."))
    }
  }
}

print("MCQRNN: Neural Network Quantile Regression")
bottom_q_NN_train = c(0.01,0.05,0.1,0.15,0.2)
bottom_q_NN_predict=seq(0.01,0.2,0.01)
top_q_NN_train=c(0.8,0.85,0.9,0.95,0.99)
top_q_NN_predict = seq(0.8,0.99,0.01)
quant_reg_mcqrnn_output <- perform_quant_reg_mcqrnn(
  orig_data,
  g_vals,
  bottom_q_NN_train = bottom_q_NN_train,
  bottom_q_NN_predict=bottom_q_NN_predict,
  top_q_NN_train=top_q_NN_train,
  top_q_NN_predict = top_q_NN_predict,
  num_hidden = 10,
  train_perc=0.9,
  user_info=FALSE
)
alternative_better_qrnn <- quant_reg_mcqrnn_output[[1]]
simp_better_qrnn <- quant_reg_mcqrnn_output[[2]]
qrnn_train_loss <- quant_reg_mcqrnn_output[[3]]
qrnn_test_loss <- quant_reg_mcqrnn_output[[4]]
qrnn_train_loss
qrnn_test_loss
for(i in 1:length(bottom_q_NN_predict)){
  print(paste0("At the alpha = ",
               bottom_q_NN_predict[i],
               " level, the alternative model is better in ",
               alternative_better_qrnn[i],
               " cases, which is a fraction of ",
               alternative_better_qrnn[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_q_NN_predict)){
  print(paste0("At the alpha = ",
               1 - top_q_NN_predict[i],
               " level, the simplified model is better in ",
               simp_better_qrnn[i],
               " cases, which is a fraction of ",
               simp_better_qrnn[i] / nrow(orig_data),
               "."))
}
if(length(bottom_q_NN_predict) == length(top_q_NN_predict)){
  for(i in 1:length(bottom_q_NN_predict)){
    if(abs(bottom_q_NN_predict[i]- (1 - top_q_NN_predict[length(top_q_NN_predict) - i+ 1])) < 0.0001){
      print(paste0("At the alpha = ",
                   bottom_q_NN_predict[i],
                   " level, the test is inconclusive in ",
                   nrow(orig_data) - (alternative_better_qrnn[i] + simp_better_qrnn[length(top_q_NN_predict) - i+ 1]),
                   " cases, which is a fraction of ",
                   (nrow(orig_data) - (alternative_better_qrnn[i] + simp_better_qrnn[length(top_q_NN_predict) - i+ 1]))/nrow(orig_data),
                   "."))
    }
  }
}

plot(bottom_quantiles,
     alternative_better_dvine_nonpar,
     type="l",
     lwd=2,
     col="darkblue",
     main="Num. Points where NCE model better",
     xlab="Quantile Level",
     ylab="Num Observations",
     ylim=range(c(alternative_better_dvine_nonpar,
                  alternative_better_dvine_par,
                  alternative_better_dvine_onepar,
                  alternative_better_lin,
                  alternative_better_qrnn)))
lines(bottom_quantiles,alternative_better_dvine_par, col="blue", lwd=2)
lines(bottom_quantiles, alternative_better_dvine_onepar, col="lightblue", lwd=2)
lines(bottom_quantiles_lin, alternative_better_lin, col="grey", lwd=2)
lines(bottom_q_NN_predict, alternative_better_qrnn, col="black", lwd=2)
# Add a legend
legend('topleft',
       legend = c('D-Vine Non-parametric', 'D-Vine Parametric', "D-Vine One Parameter", 'Linear QR', "MCQRNN"),
       col = c('darkblue', 'blue', 'lightblue', "grey", "black"),
       lty = 1,           # Line type (1 for solid line)
       lwd = 2)

plot(top_quantiles,
     simp_better_dvine_nonpar,
     type="l",
     main="Num. Points where simplified model better",
     xlab="Quantile Level",
     ylab="Num Observations",
     ylim=range(c(simp_better_dvine_nonpar,
                  simp_better_dvine_par,
                  simp_better_dvine_onepar,
                  simp_better_lin,
                  simp_better_qrnn)),
     lwd=2)
lines(top_quantiles,simp_better_dvine_par, col="blue",lwd=2)
lines(top_quantiles,simp_better_dvine_onepar, col="lightblue",lwd=2)
lines(top_quantiles_lin, simp_better_lin, col="grey",lwd=2)
lines(top_q_NN_predict, simp_better_qrnn, col="black",lwd=2)
# Add a legend
legend('topleft',
       legend = c('D-Vine Non-parametric', 'D-Vine Parametric', "D-Vine One Parameter",'Linear QR', "MCQRNN"),
       col = c('darkblue', 'blue', 'lightblue', "grey", "black"),
       lty = 1,           # Line type (1 for solid line)
       lwd = 2)

# Using ggplot2:
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

# Combine data into a long-format tibble
plot_df <- tibble(
  quantile_level = c(bottom_quantiles, bottom_quantiles, bottom_quantiles, bottom_quantiles_lin, bottom_q_NN_predict),
  observations = c(
    alternative_better_dvine_nonpar,
    alternative_better_dvine_par,
    alternative_better_dvine_onepar,
    alternative_better_lin,
    alternative_better_qrnn
  ),
  model = factor(rep(c("D-Vine Non-parametric", "D-Vine Parametric", "D-Vine One-parametric",
                       "Linear QR", "MCQRNN"),
                     times = c(length(bottom_quantiles), length(bottom_quantiles),
                               length(bottom_quantiles), length(bottom_quantiles_lin),
                               length(bottom_q_NN_predict))))
)
bottom_quantiles_plot <- ggplot(plot_df, aes(x = quantile_level, y = observations, color = model)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Num. Points where NCE model better",
    x = "Quantile Level",
    y = "Num Observations",
    color = "Model"
  ) +
  scale_color_manual(values = c(
    "D-Vine Non-parametric" = "darkblue",
    "D-Vine Parametric" = "blue",
    "D-Vine One-parametric" = "lightblue",
    "Linear QR" = "grey",
    "MCQRNN" = "black"
  )) +
  theme_minimal(base_size = 14)

bottom_quantiles_plot
ggsave(
  filename = "BottomQuantilesLinearExample20250625.png",
  plot = bottom_quantiles_plot,
  width = 10,
  height = 6,
  dpi = 300,
  bg="white"
)

# For top quantiles:
# Combine data into a long-format tibble
plot_df_top <- tibble(
  quantile_level = c(top_quantiles, top_quantiles, top_quantiles, top_quantiles_lin, top_q_NN_predict),
  observations = c(
    simp_better_dvine_nonpar,
    simp_better_dvine_par,
    simp_better_dvine_onepar,
    simp_better_lin,
    simp_better_qrnn
  ),
  model = factor(rep(c("D-Vine Non-parametric", "D-Vine Parametric", "D-Vine One-parametric",
                       "Linear QR", "MCQRNN"),
                     times = c(length(top_quantiles), length(top_quantiles),
                               length(top_quantiles), length(top_quantiles_lin),
                               length(top_q_NN_predict))))
)
top_quantiles_plot <- ggplot(plot_df_top, aes(x = quantile_level, y = observations, color = model)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Num. Points where simplified model better",
    x = "Quantile Level",
    y = "Num Observations",
    color = "Model"
  ) +
  scale_color_manual(values = c(
    "D-Vine Non-parametric" = "darkblue",
    "D-Vine Parametric" = "blue",
    "D-Vine One-parametric" = "lightblue",
    "Linear QR" = "grey",
    "MCQRNN" = "black"
  )) +
  theme_minimal(base_size = 14)

top_quantiles_plot
ggsave(
  filename = "TopQuantilesLinearExample20250625.png",
  plot = top_quantiles_plot,
  width = 10,
  height = 6,
  dpi = 300,
  bg="white"
)

# save g_vals
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
csv_file_name <- paste0("data/g_values_",ncol(orig_data),"d_", current_date, ".csv")
# Save as CSV
write.csv(g_vals, file = csv_file_name, row.names = FALSE)
print(paste("Data saved to", csv_file_name))



# Test true density
tau_upper = 0.92
tau_lower=-0.92
struct_mat <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
family_test <- list(list("frank", "clayton","gaussian","frank"),
                    list("frank","gaussian","joe"),
                    list("gaussian", "gumbel"),
                    list("gaussian"))
params_test <- list(c(ktau_to_par(family=family_test[[1]][[1]], tau=-0.2)),
                    c(ktau_to_par(family=family_test[[1]][[2]], tau=0.3)),
                    c(ktau_to_par(family=family_test[[1]][[3]], tau=-0.1)),
                    c(ktau_to_par(family=family_test[[1]][[4]], tau=0.1)))
param_cond_funcs_test <- list(
  list(u_to_param_linear(c(1),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper),
       u_to_param_linear(c(1),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper),
       u_to_param_linear(c(1),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper)),
  list(u_to_param_linear(c(0.7,0.3),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper),
       u_to_param_linear(c(0.4,0.6),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper)),
  list(u_to_param_linear(c(0.2,0.5,0.3),
                         tau_lower=tau_lower,
                         tau_upper=tau_upper)))
true_log_lik <- log_likelihood_non_simplified(u_data =orig_data,
                                          struct = struct_mat,
                                          families=family_test,
                                          params = params_test,
                                          param_cond_funcs = param_cond_funcs_test,
                                          rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)),
                                          return_vector=TRUE)

