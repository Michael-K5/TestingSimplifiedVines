# Fit a simplified vine copula and train a classifier on the data, to distinguish between
# Non-Simplified and Simplified Vine Copula data.
library(rvinecopulib)
library(keras)
library(tensorflow)
source("R/classifier_methods.R")
source("R/plotting_methods.R")
# Parameters to determine, which data to load.
last_data_simulation_date <- "2025-06-25"
data_dim <- 5
# fraction of noise to true samples (to determine how many noise samples to create)
nu <- 1
# load data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"d_",last_data_simulation_date,".csv")
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames
num_rows <- nrow(orig_data)
#pairs_copula_data(orig_data)

# fit a simplified vine
fitted_vine <-vinecop(orig_data,family_set="parametric")
print.data.frame(summary(fitted_vine),digit=2)
fitted_vine_details_df <- summary(fitted_vine)
fitted_vine_details_df
readr::write_csv(fitted_vine_details_df,
          file="simplified_vine_20250625_details.csv")
knitr::kable(fitted_vine_details_df, format="latex", booktabs="TRUE")
# save the fitted_vine
current_date <- Sys.Date()
copula_path <- paste0("models/copula_",data_dim,"d_", current_date,".rds")
saveRDS(fitted_vine, file = copula_path)
# simulate from the simplified vine, to train a classifier later.
num_samples <- num_rows * nu# number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
# plot at most first 10 k for faster rendering
#pairs_copula_data(simplified_samples[(1:min(10000, nrow(num_samples))),])
simplified_pairs_plot <- copula_pairs_ggplot(simplified_samples[(1:min(10000, nrow(num_samples))),])
simplified_pairs_plot <- ggplotify::as.ggplot(simplified_pairs_plot)
simplified_pairs_plot
# ggplot2::ggsave(
#   filename = "simplified_samples_5d_20250625_plot.png",
#   plot = simplified_pairs_plot,
#   width = 6,
#   height = 6,
#   dpi = 300,
#   bg="white"
# )
split_output <- train_test_split(orig_data, simplified_samples, stratified=TRUE)
x_train <- split_output[[1]]
x_test <- split_output[[2]]
y_train <- split_output[[3]]
y_test <- split_output[[4]]

#nrow(x_train[y_train==1,])
#nrow(x_test[y_test==1,])
#nrow(x_train)
model <- build_model(
    input_dim=5,
    hidden_units=c(32,16), # 2 hidden layers with 20 and 10 units respectively.
    initial_lr = 0.01,
    use_tanh=FALSE, # Use leaky_relu, not tanh
    leaky_relu_alpha=0.1)
start_time <- Sys.time()
train_model_output <- train_model(
  model=model,
  x_train=x_train,
  y_train=y_train,
  lr_schedule=lr_schedule_fun,
  num_epochs=200
)
end_time <- Sys.time()
print(end_time - start_time)
model <- train_model_output[[1]]
history <- train_model_output[[2]]

# training and validation loss and accuracy, learning rate
# plot(history)
# remove learning rate from plot
plot(history, metrics = c("loss", "binary_accuracy_from_logits"))
# evaluate the model on the test set
loss_and_metrics <- model %>% evaluate(x_test, y_test)
print(loss_and_metrics[["loss"]])
print(loss_and_metrics[["binary_accuracy_from_logits"]])
loss_and_metrics_train <- model %>% evaluate(x_train, y_train)
print(loss_and_metrics_train[["loss"]])
print(loss_and_metrics_train[["binary_accuracy_from_logits"]])
print(paste0(
  "Base Accuracy for Prior: ",
  nrow(simplified_samples) / (nrow(simplified_samples) + nrow(orig_data))
  )
)
int_val <- compute_integral(model, fitted_vine, n_samples=20000, nu=nu,data_dim_if_unif=ncol(orig_data),user_info=TRUE)
int_val
# save the model for reusing it later
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
model_file_name <- paste0("models/NN_", ncol(x_train), "d_",current_date, ".keras")
#save_model_hdf5(model, filepath = model_file_name)
keras$Model$save(model, filepath=model_file_name)

log_lik_simp_train_orig <- sum(log(dvinecop(x_train[y_train==1,], fitted_vine)))
log_lik_NN_train_orig <- sum(log(non_param_cop(
  model=model,
  fitted_vine=fitted_vine,
  obs=x_train[y_train==1,],
  nu=nu)))
log_lik_simp_test_orig <- sum(log(dvinecop(x_test[y_test==1,], fitted_vine)))
log_lik_NN_test_orig <- sum(log(non_param_cop(
  model=model,
  fitted_vine=fitted_vine,
  obs=x_test[y_test==1,],
  nu=nu)))
NN_num_params <- count_NN_params(weights=model$weights)
simp_cop_num_params <- get_num_cop_params(fitted_vine)
deg_free <- NN_num_params - simp_cop_num_params
LRT_stat_train <- -2*(as.numeric(log_lik_simp_train_orig) - as.numeric(log_lik_NN_train_orig))
LRT_stat_test <- -2*(as.numeric(log_lik_simp_test_orig) - as.numeric(log_lik_NN_test_orig))
# compute p(chisq(df=deg_free) > LRT_stat) (for train and test set)
p_value_train <- pchisq(LRT_stat_train, df = deg_free, lower.tail = FALSE)
p_value_test <- pchisq(LRT_stat_test, df = deg_free, lower.tail = FALSE)
AIC_NN <- 2*(NN_num_params+simp_cop_num_params) - 2* (log_lik_NN_train_orig + log_lik_NN_test_orig)
BIC_NN <- log(nrow(orig_data))*(NN_num_params+ simp_cop_num_params) - 2* (log_lik_NN_train_orig + log_lik_NN_test_orig)
AIC_simp <- 2*simp_cop_num_params - 2*(log_lik_simp_train_orig + log_lik_simp_test_orig)
BIC_simp <- log(nrow(orig_data))*simp_cop_num_params - 2*(log_lik_simp_train_orig + log_lik_simp_test_orig)
eval_data_frame <- data.frame(
  MC_Integral = int_val,
  train_accuracy= round(loss_and_metrics_train[["binary_accuracy_from_logits"]],4),
  train_loss = round(loss_and_metrics_train[["loss"]],4),
  test_accuracy = round(loss_and_metrics[["binary_accuracy_from_logits"]],4),
  test_loss = round(loss_and_metrics[["loss"]],4),
  log_lik_simplified_train = round(log_lik_simp_train_orig,2),
  log_lik_NN_train = round(log_lik_NN_train_orig,2),
  log_lik_simplified_test = round(log_lik_simp_test_orig,2),
  log_lik_NN_test = round(log_lik_NN_test_orig,2),
  AIC_NN = AIC_NN,
  BIC_NN=BIC_NN,
  AIC_simp = AIC_simp,
  BIC_simp = BIC_simp,
  NN_num_params = NN_num_params,
  simp_cop_num_params = simp_cop_num_params,
  LRT_stat_train = LRT_stat_train,
  LRT_stat_test = LRT_stat_test,
  p_value_train = p_value_train,
  p_value_test = p_value_test
)
readr::write_csv(eval_data_frame, "ExampleRunEvaluationMetrics20250625.csv")
print(p_value_train)
print(p_value_test)
print(AIC_NN)
print(AIC_simp)
print(BIC_NN)
print(BIC_simp)
print(log_lik_NN_train_orig)
print(log_lik_simp_train_orig)
print(log_lik_NN_test_orig)
print(log_lik_simp_test_orig)

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
