# Fit a simplified vine copula and train a classifier on the data, to distinguish between
# Non-Simplified and Simplified Vine Copula data.
library(rvinecopulib)
library(keras)
library(tensorflow)
# Parameters to determine, which data to load.
last_data_simulation_date <- "2025-05-08"
data_dim <- "5"
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
num_samples <- num_rows # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
#pairs_copula_data(simplified_samples)

# If required: have same ratio in train and test sets
train_test_split <- function(orig_data, simplified_data, train_perc=0.8, stratified=FALSE){
  # cast labels as matrix (to ensure compatibility with tensorflow)
  labels <- as.matrix(rep(1L, nrow(orig_data)))
  simplified_labels <- as.matrix(rep(0L, nrow(simplified_data)))
  if(stratified){
    # extract the training and testing samples from the original data
    num_train_orig <- floor(train_perc * nrow(orig_data))
    train_orig_idx <- sample(nrow(orig_data),num_train_orig)
    x_train_orig <- orig_data[train_orig_idx,]
    x_test_orig <- orig_data[-train_orig_idx,]
    y_train_orig <- as.matrix(labels[train_orig_idx,])
    y_test_orig <- as.matrix(labels[-train_orig_idx,])
    # extract the training and testing samples from the simplified data
    num_train_simp <- floor(train_perc*nrow(simplified_samples))
    train_simp_idx <- sample(nrow(simplified_samples), num_train_simp)
    x_train_simp <- simplified_samples[train_simp_idx,]
    y_train_simp <- as.matrix(simplified_labels[train_simp_idx,])
    x_test_simp <- simplified_samples[-train_simp_idx,]
    y_test_simp <- as.matrix(simplified_labels[-train_simp_idx,])
    # combine the data into one dataset and shuffle
    x_train <- rbind(x_train_orig, x_train_simp)
    y_train <- rbind(y_train_orig, y_train_simp)
    shuffle_idx_train <- sample(nrow(x_train), nrow(x_train))
    x_train <- x_train[shuffle_idx_train,]
    y_train <- y_train[shuffle_idx_train,]
    x_test <- rbind(x_test_orig, x_test_simp)
    y_test <- rbind(y_test_orig, y_test_simp)
    shuffle_idx_test <- sample(nrow(x_test), nrow(x_test))
    x_test <- x_test[shuffle_idx_test,]
    y_test <- y_test[shuffle_idx_test,]
    return(list(x_train,
                x_test,
                y_train,
                y_test))
  } else {
    # combine data
    classifier_data <- rbind(orig_data, simplified_samples)
    classifier_labels <- rbind(labels, simplified_labels)
    # perform train test split
    num_train <- floor(train_perc * nrow(classifier_data))
    train_idx <- sample(nrow(classifier_data),num_train)
    x_train <- classifier_data[train_idx,]
    y_train <- as.matrix(classifier_labels[train_idx,])
    x_test <- classifier_data[-train_idx,]
    y_test <- as.matrix(classifier_labels[-train_idx,])
    return(list(x_train,
                x_test,
                y_train,
                y_test))
  }
}

split_output <- train_test_split(orig_data, simplified_samples)
x_train <- split_output[[1]]
x_test <- split_output[[2]]
y_train <- split_output[[3]]
y_test <- split_output[[4]]

# definition of the model
model <- keras_model_sequential() %>%
  #layer_dense(units = 20, activation = 'relu', input_shape = ncol(x_train)) %>%
  layer_dense(units = 32, input_shape = ncol(x_train)) %>%
  layer_activation_leaky_relu(alpha = 0.1) %>%
  #layer_dropout(rate = 0.3) %>% # Dropout layer (for regularization, if necessary)
  #layer_dense(units = 10, activation = 'relu') %>%
  layer_dense(units = 16) %>%
  layer_activation_leaky_relu(alpha = 0.1) %>%
  layer_dense(units = 1, activation = 'sigmoid') # sigmoid activation, for binary classification

# learning rate scheduler: halves the learning rate every 30 epochs
lr_schedule <- function(epoch, lr) {
  if((epoch + 1) %% 10 == 0){
    return(lr/2)
  } else {
    return(lr)
  }
}

# create the keras object required for learning rate scheduling
lr_scheduler <- callback_learning_rate_scheduler(schedule = lr_schedule)

# compile the model, define optimizer, loss and metrics
model %>% compile(
  optimizer = keras$optimizers$Adam(learning_rate=1e-2),
  loss = keras$losses$BinaryCrossentropy(),
  metrics = keras$metrics$BinaryAccuracy()
)
# show model summary
model %>% summary

# train the model
history <- model %>% fit(
  x_train, y_train,
  epochs = 50,
  batch_size = 128,
  validation_split=0.2, # use 20 percent of training data as validation data
  verbose = 1, # 0 for slightly faster training (no output), 1 to observe progress while training
  callbacks=list(lr_scheduler)
)
plot(history)

# evaluate the model on the test set
loss_and_metrics <- model %>% evaluate(x_test, y_test)
print(loss_and_metrics)

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
