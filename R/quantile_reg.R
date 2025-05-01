# Script for performing a quantile regression to test for the simmplifying assumption.
# Also contains a function for evaluating a non parametric copula estimate

# Parameters to determine, which model and data to load.
library(vinereg)
library(keras)
library(rvinecopulib)
last_data_simulation_date <- "2025-04-29"
last_train_date <- "2025-04-29"
data_dim <- "4"
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


# compute the value of the non-parametric copula obtained from the simplified fit,
# together with the classifier
non_param_cop <- function(obs){
  predictions <- model %>% predict(obs)
  r_vals <- predictions/(1-predictions)
  c_np <- exp(r_vals) * dvinecop(obs, fitted_cop)
  return(c_np)
}
temp <- non_param_cop(orig_data)
head(temp)

# make predictions on the observed Data (orig_data)
head(orig_data)
predictions <- model %>% predict(orig_data)
# how many samples the model recognizes as non-simplified
length(predictions[predictions> 0.5]) / length(predictions)
# for binary predictions: binary_predictions <- ifelse(predictions > 0.5, 1, 0)


# here a classifier p(u) was trained on the data. to get the values r(u),
# with r(u)/(1+r(u)) = p(u), we need to compute r(u) = p(u)/(1-p(u))
r_vals <- predictions / (1-predictions)
head(r_vals)
head(predictions)

# quantile regression
obs_data <- data.frame(orig_data)
qreg_data <- cbind(r_vals, obs_data)
q_reg <- vinereg(r_vals ~ . ,family_set="parametric", data=qreg_data)
quantiles <- predict(q_reg, newdata=obs_data, alpha=c(0.10,0.90))
alternative_better <- 0
simp_better <- 0
for( i in 1:nrow(quantiles)){
  if(quantiles[i,1] > 1){
    alternative_better <- alternative_better + 1
  } else if (quantiles[i,2] < 1){
    simp_better <- simp_better +1
  }
}
alternative_better
simp_better
nrow(quantiles) - (alternative_better + simp_better)
alternative_better / nrow(quantiles)
# TODO: Test and Format output better.
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
csv_file_name <- paste0("data/r_values_",ncol(orig_data),"d_", current_date, ".csv")
# Save as CSV
write.csv(r_vals, file = csv_file_name, row.names = FALSE)
print(paste("Data saved to", csv_file_name))
