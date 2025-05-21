# Script for performing a quantile regression to test for the simmplifying assumption.
# Also contains a function for evaluating a non parametric copula estimate
#library(vinereg)
library(keras)
library(rvinecopulib)
# Parameters to determine, which model and data to load.
last_data_simulation_date <- "2025-05-08"
last_train_date <- "2025-05-08"
data_dim <- "5"
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


#' compute the value of the non-parametric copula obtained from the simplified fit,
#' together with the classifier
#' @param model: The classifier trained to distinguish non-simplified data
#' with labels 1 from simplified data with labels 0
#' @param fitted_vine: The vine copula fitted to the observed data.
#' @param obs: the observations, for which the non-parametric copula should be created.
#' @param nu: The fraction of noise to real samples
#' @returns The non-parametric copula evaluated at the points obs.
non_param_cop <- function(model, fitted_vine, obs, nu=1){
  predictions <- model %>% predict(obs)
  c_np <- nu * predictions/(1-predictions) * dvinecop(obs, fitted_vine)
  return(c_np)
}

#' monte carlo integral of non_param_cop (using importance sampling)
#' @param model: The classifier trained to distinguish non-simplified data
#' with labels 1 from simplified data with labels 0
#' @param fitted_vine: The vine copula fitted to the observed data.
#' @param n_samples: Number of samples to create to compute the integral
#' @param user_info: Whether to display an information, which steps are currently running.
#' @returns The monte carlo integral approximation.
compute_integral <- function(model, fitted_vine, n_samples, user_info=FALSE){
  if (user_info){
    print("Start noise sampling")
  }
  samples <- rvinecop(n_samples, fitted_vine)
  if(user_info){
    print("Noise samples created")
    print("Evaluating noise density")
  }
  p_simp <- dvinecop(samples, fitted_cop)
  if(user_info){
    print("evaluating neural network output")
  }
  p_non_param <- non_param_cop(fitted_vine=fitted_vine, model=model, obs=samples)
  return(mean(p_non_param/p_simp))
}
# temp <- non_param_cop(orig_data)
# int_val <- compute_integral(100000)
# int_val
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
obs_data <- data.frame(orig_data)
qreg_data <- cbind(g_vals, obs_data)
q_reg <- vinereg(g_vals ~ . ,family_set="parametric", data=qreg_data)
quantiles <- predict(q_reg, newdata=obs_data, alpha=c(0.1,0.9))
alternative_better <- 0
simp_better <- 0
for( i in 1:nrow(quantiles)){
  if(quantiles[i,1] > 0){
    alternative_better <- alternative_better + 1
  } else if (quantiles[i,2] < 0){
    simp_better <- simp_better +1
  }
}
alternative_better
simp_better
nrow(quantiles) - (alternative_better + simp_better)
alternative_better / nrow(quantiles)
# head(quantiles)
# temp <- 0
# for (i in 1:nrow(quantiles)){
#   if (quantiles[i,1] < r_vals[i]){
#     temp <- temp +1
#   }
# }
# temp


# save r_vals
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
csv_file_name <- paste0("data/g_values_",ncol(orig_data),"d_", current_date, ".csv")
# Save as CSV
write.csv(g_vals, file = csv_file_name, row.names = FALSE)
print(paste("Data saved to", csv_file_name))
