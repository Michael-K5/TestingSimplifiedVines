# Fit a simplified vine copula and train a classifier on the data, to distinguish between
# Non-Simplified and Simplified Vine Copula data.
library(rvinecopulib)
library(keras)
library(tensorflow)
source("R/classifier_methods.R")
# Parameters to determine, which data to load.
last_data_simulation_date <- "2025-06-02"
data_dim <- 5
# fraction of noise to true samples (to determine how many noise samples to create)
nu <- 4
# load data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"d_",last_data_simulation_date,".csv")
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames
num_rows <- nrow(orig_data)
#pairs_copula_data(orig_data)

# fit a simplified vine
fitted_vine <-vinecop(orig_data,family_set="parametric")
print.data.frame(summary(fitted_vine),digit=2)
# save the fitted_vine
current_date <- Sys.Date()
copula_path <- paste0("models/copula_",data_dim,"d_", current_date,".rds")
saveRDS(fitted_vine, file = copula_path)
# simulate from the simplified vine, to train a classifier later.
num_samples <- num_rows * nu# number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
# plot at most first 10 k for faster rendering
pairs_copula_data(simplified_samples[1:min(10000, num_samples),])

split_output <- train_test_split(orig_data, simplified_samples, stratified=TRUE)
x_train <- split_output[[1]]
x_test <- split_output[[2]]
y_train <- split_output[[3]]
y_test <- split_output[[4]]

model <- build_model(
    input_dim=5,
    hidden_units=c(20,10), # 2 hidden layers with 20 and 10 units respectively.
    initial_lr = 0.01,
    use_tanh=FALSE, # Use leaky_relu, not tanh
    leaky_relu_alpha=0.1)
train_model_output <- train_model(
  model=model,
  x_train=x_train,
  y_train=y_train,
  lr_schedule=lr_schedule_fun,
  num_epochs=200
)
model <- train_model_output[[1]]
history <- train_model_output[[2]]

plot(history) # training and validation loss and accuracy, learning rate
# remove learning rate from plot
plot(history, metrics = c("loss", "binary_accuracy"))
# evaluate the model on the test set
loss_and_metrics <- model %>% evaluate(x_test, y_test)
print(loss_and_metrics)
print(paste0(
  "Base Accuracy for Prior: ",
  nrow(simplified_samples) / (nrow(simplified_samples) + nrow(orig_data))
  )
)
int_val <- compute_integral(model, fitted_vine, n_samples=50000, nu=4,data_dim_if_unif=5,user_info=TRUE)
int_val
# save the model for reusing it later
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
model_file_name <- paste0("models/NN_", ncol(x_train), "d_",current_date, ".keras")
#save_model_hdf5(model, filepath = model_file_name)
keras$Model$save(model, filepath=model_file_name)
# The syntax for loading and using this model is:
# loaded_model <- load_model_hdf5("models/NN_5d_2025-04-29.keras")
# # example for predictions, further training or evaluating the loaded model.
# predictions <- predict(loaded_model, orig_data)
# history <- loaded_model %>% fit(
#   x_train, y_train,
#   epochs = 2,
#   batch_size = 100,
#   validation_split=0.2, # use 20 percent of training data as validation data
#   verbose = 1 # 0 for slightly faster training (no output), 1 to observe progress while training
# )
# loaded_model %>% evaluate(x_test, y_test)
