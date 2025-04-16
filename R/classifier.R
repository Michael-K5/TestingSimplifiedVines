# Fit a simplified vine copula and train a classifier on the data
library(rvinecopulib)
library(keras)
library(tensorflow)
# load the original data
csv_filename <- "data/non_simplified_sim_2025-04-16.csv"
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames
num_rows <- nrow(orig_data)
# cast labels as matrix (otherwise this does not work with tensorflow)
labels <- as.matrix(rep(1L, num_rows))
#pairs_copula_data(orig_data)

# fit a simplified vine
fitted_vine <-vinecop(orig_data,family_set="onepar")
print.data.frame(summary(fitted_vine),digit=2)
num_samples <- num_rows # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
#pairs_copula_data(simplified_samples)
simplified_labels <- as.matrix(rep(0L, num_samples))

classifier_data <- rbind(orig_data, simplified_samples)
classifier_labels <- rbind(labels, simplified_labels)

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
  epochs = 500,
  batch_size = 100,
  validation_split=0.2, # use 20 percent of training data as validation data
  verbose = 1 # 0 for slightly faster training (no output), 1 to observe progress while training
)
plot(history)

# evaluate the model on the test set
loss_and_metrics <- model %>% evaluate(x_test, y_test)

# compute the value of the non-parametric copula obtained from the simplified fit,
# together with the classifier
non_param_cop <- function(obs){
  predictions <- model %>% predict(obs)
  r_vals <- predictions/(1-predictions)
  c_np <- exp(r_vals) * dvinecop(obs, fitted_vine)
  return(c_np)
}

# make predictions on the observed Data (orig_data)
head(orig_data)
predictions <- model %>% predict(orig_data)
head(predictions)
# for binary predictions: binary_predictions <- ifelse(predictions > 0.5, 1, 0)

# here a classifier p(u) was trained on the data. to get the values r(u),
# with r(u)/(1+r(u)) = p(u), we need to compute r(u) = p(u)/(1-p(u))
r_vals <- predictions / (1-predictions)

head(r_vals)

# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()

# Construct the file name with the date
csv_file_name <- paste0("data/r_values_", current_date, ".csv")

# Save as CSV
write.csv(r_vals, file = csv_file_name, row.names = FALSE)

print(paste("Data saved to", csv_file_name, "\n"))

