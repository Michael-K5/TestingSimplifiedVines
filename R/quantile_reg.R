# Script for performing a quantile regression to test for the simmplifying assumption.
# Also contains a function for evaluating a non parametric copula estimate
#library(vinereg)
source("R/classifier_methods.R")
library(keras)
library(rvinecopulib)
# Parameters to determine, which model and data to load.
last_data_simulation_date <- "2025-05-23"
last_train_date <- "2025-05-23"
data_dim <- "5"
nu <- 5 # T_n/T_c, the ratio of noise to observed samples during training
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
plot(cor_facs)
remove_top <- floor(0.01 * length(cor_facs))
cor_facs_no_outliers <- (sort(cor_facs)[1:(length(cor_facs) - remove_top)])
plot(cor_facs_no_outliers)
hist(cor_facs_no_outliers)
sum(cor_facs_no_outliers < 1)
sum(cor_facs_no_outliers > 2)

# temp <- non_param_cop(orig_data)
int_val <- compute_integral(model=model,
                            fitted_vine = fitted_cop,
                            nu=nu,
                            data_dim_if_unif = as.integer(data_dim),
                            n_samples=10000,
                            user_info=TRUE)
int_val
# temp_norm <- temp / int_val
#
# make predictions on the observed Data (orig_data)
head(orig_data)
predictions <- model %>% predict(orig_data)
# how many samples the model recognizes as non-simplified
length(predictions[predictions> 0.5]) / length(predictions)
# for binary predictions: binary_predictions <- ifelse(predictions > 0.5, 1, 0)


# here a classifier p(u) was trained on the data. to get the values g(u),
# as defined in the thesis, g(u) = log(\nu \cdot p(u)/(1-p(u))) needs to be computed,
# where nu = T_n/T_c, the fraction of number of noise samples to observed samples
g_vals <- log(nu * predictions / (1-predictions)) # G(u, \eta) from the thesis

# quantile regression:
# Test if the 10 percent quantile of G is >0, then the alternative model is considered better
# If the 90 percent quantile of G is < 0, then the simplified model is better (unlikely)
bottom_q_levs <- c(0.05, 0.1)
top_q_levs <- c(0.9,0.95)
output_qreg <- perform_quant_reg(
    g_vals,
    nu=5,
    orig_data,
    family_set_name = "parametric",
    bottom_quantile_levels = bottom_q_levs,
    top_quantile_levels = top_q_levs)

alternative_better <- output_qreg[[1]]
simp_better <- output_qreg[[2]]
for(i in 1:length(bottom_q_levs)){
  print(paste0("At the alpha = ",
               bottom_q_levs[i],
               " level, the alternative model is better in ",
               alternative_better[i],
               " cases, which is a fraction of ",
               alternative_better[i] / nrow(orig_data),
               "."))
}
for(i in 1:length(top_q_levs)){
  print(paste0("At the alpha = ",
               1 - top_q_levs[i],
               " level, the simplified model is better in ",
               simp_better[i],
               " cases, which is a fraction of ",
               simp_better[i] / nrow(orig_data),
               "."))
}
if(length(bottom_q_levs) == length(top_q_levs)){
  for(i in 1:length(bottom_q_levs)){
    if(abs(bottom_q_levs[i]- (1 - top_q_levs[length(top_q_levs) - i+ 1])) < 0.0001){
      print(paste0("At the alpha = ",
                   bottom_q_levs[i],
                   " level, the test is inconclusive in ",
                   nrow(orig_data) - (alternative_better[i] + simp_better[length(top_q_levs) - i+ 1]),
                   " cases, which is a fraction of ",
                   (nrow(orig_data) - (alternative_better[i] + simp_better[length(top_q_levs) - i+ 1]))/nrow(orig_data),
                   "."))
    }
  }
}
# head(quantiles)
# temp <- 0
# for (i in 1:nrow(quantiles)){
#   if (quantiles[i,1] < r_vals[i]){
#     temp <- temp +1
#   }
# }
# temp


# save g_vals
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
csv_file_name <- paste0("data/g_values_",ncol(orig_data),"d_", current_date, ".csv")
# Save as CSV
write.csv(g_vals, file = csv_file_name, row.names = FALSE)
print(paste("Data saved to", csv_file_name))
