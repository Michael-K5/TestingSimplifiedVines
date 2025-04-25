# Fit a simplified vine copula and train a classifier on the data, to distinguish between
# Non-Simplified and Simplified Vine Copula data.
# TODO: Implement a cross validation approach, instead of just one test set.
library(rvinecopulib)
library(keras)
library(tensorflow)
# Parameters to determine, which data to load.
last_data_simulation_date <- "2025-04-25"
data_dim <- "5"
# load data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"d_",last_data_simulation_date,".csv")
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames
num_rows <- nrow(orig_data)
# cast labels as matrix (otherwise this does not work with tensorflow)
labels <- as.matrix(rep(1L, num_rows))
#pairs_copula_data(orig_data)

# fit a simplified vine
fitted_vine <-vinecop(orig_data,family_set="onepar")
print.data.frame(summary(fitted_vine),digit=2)
# save the fitted_vine
current_date <- Sys.Date()
copula_path <- paste0("models/copula_",data_dim,"d_", current_date,".rds")
saveRDS(fitted_vine, file = copula_path)
# simulate from the simplified vine, to train a classifier later.
num_samples <- num_rows # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
#pairs_copula_data(simplified_samples)
simplified_labels <- as.matrix(rep(0L, num_samples))

classifier_data <- rbind(orig_data, simplified_samples)
classifier_labels <- rbind(labels, simplified_labels)
# TODO: make the following a k-fold cross validation instead (k=5)
# perform train test split
train_perc <- 0.8
num_train <- floor(train_perc * nrow(classifier_data))
train_idx <- sample(nrow(classifier_data),num_train)
x_train <- classifier_data[train_idx,]
y_train <- as.matrix(classifier_labels[train_idx,])
x_test <- classifier_data[-train_idx,]
y_test <- as.matrix(classifier_labels[-train_idx,])
# head(x_train)
# head(x_test)
# head(y_train)
# head(y_test)

# definition of the model
model <- keras_model_sequential() %>%
  layer_dense(units = 20, activation = 'relu', input_shape = ncol(x_train)) %>%
  #layer_dropout(rate = 0.3) %>% # Dropout layer (for regularization, if necessary)
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dense(units = 1, activation = 'sigmoid') # sigmoid activation, for binary classification


# compile the model, define optimizer, loss and metrics
model %>% compile(
  optimizer = 'adam',
  loss = 'binary_crossentropy', #
  metrics = c('accuracy')
)
# show model summary
model %>% summary

# train the model
history <- model %>% fit(
  x_train, y_train,
  epochs = 300,
  batch_size = 100,
  validation_split=0.2, # use 20 percent of training data as validation data
  verbose = 1 # 0 for slightly faster training (no output), 1 to observe progress while training
)
plot(history)

# evaluate the model on the test set
loss_and_metrics <- model %>% evaluate(x_test, y_test)

# save the model for reusing it later
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
model_file_name <- paste0("models/NN_", ncol(x_train), "d_",current_date, ".keras")
save_model_hdf5(model, filepath = model_file_name)

# The syntax for loading and using this model is:
# loaded_model <- load_model_hdf5("models/NN_3d_2025-04-16.keras")
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
