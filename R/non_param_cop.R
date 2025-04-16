# Script for evaluating a non parametric copula estimate
# Parameters to determine, which model and data to load.
last_data_simulation_date <- "2025-04-16"
last_train_date <- "2025-04-16"
data_dim <- "3d"
# load the Neural Network Classifier
model_path <- paste0("models/NN_", data_dim, "_", last_train_date, ".keras")
model <- load_model_hdf5(model_path)
# load the fitted simplified vine
cop_path <- paste0("models/copula_",data_dim,"_", last_train_date,".rds")
fitted_cop <- readRDS(file = cop_path)
# load the original data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"_",last_data_simulation_date,".csv")
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
head(predictions)
# for binary predictions: binary_predictions <- ifelse(predictions > 0.5, 1, 0)

# here a classifier p(u) was trained on the data. to get the values r(u),
# with r(u)/(1+r(u)) = p(u), we need to compute r(u) = p(u)/(1-p(u))
r_vals <- predictions / (1-predictions)
head(r_vals)
# TODO: perform distribution regression analysis on the r values
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()

# Construct the file name with the date
csv_file_name <- paste0("data/r_values_", current_date, ".csv")

# Save as CSV
write.csv(r_vals, file = csv_file_name, row.names = FALSE)

print(paste("Data saved to", csv_file_name, "\n"))



