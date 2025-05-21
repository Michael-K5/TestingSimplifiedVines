library(rvinecopulib) # copulas
library(keras) # Machine Learning and Neural Networks
library(tensorflow) # Neural Networks

#' Performs a train test split
#' Gives all datapoints in orig_data a label of 1, all datapoints in simplified data a label of 0
#' @param orig_data: The real dataset, which gets labels of 1
#' @param simplified_data: Synthetic dataset, simulated from a simplified
#' vine copula, gets labels of 0
#' @param train_perc: Number between 0 and 1, defaults to 0.8.
#' Determins what fraction of the data is used for training.
#' @param stratified: TRUE or FALSE, defaults to FALSE
#' Whether the share of points from orig data should be exactly the same in train and test set.
#' @returns A list with 4 entries, x_train, x_test, y_train, y_test in that order.
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
    num_train_simp <- floor(train_perc*nrow(simplified_data))
    train_simp_idx <- sample(nrow(simplified_data), num_train_simp)
    x_train_simp <- simplified_data[train_simp_idx,]
    y_train_simp <- as.matrix(simplified_labels[train_simp_idx,])
    x_test_simp <- simplified_data[-train_simp_idx,]
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
    classifier_data <- rbind(orig_data, simplified_data)
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

#' learning rate scheduler: halves the learning rate every 30 epochs
#' argument to a keras function
#' @param epoch: current epoch
#' @param lr: current learning rate
#' @returns lr: the updated learning rate
lr_schedule_fun <- function(epoch, lr) {
  if((epoch + 1) %% 30 == 0){
    return(lr/2)
  } else {
    return(lr)
  }
}
#' Defines and compiles a neural network for binary classification.
#' Activation function is leaky_relu, output is one layer with sigmoid activation
#' @param input_dim: dimension of the input, defaults to 5
#' @param hidden_units: vector containing the number of units per layer.
#' Defaults to c(20,10)
#' @param initial_lr: initial learing rate to use, defaults to 0.01
#' (íf no lr_scheduler is used during training, this stays the same during
#' the whole training process)
#' @param leaky_relu_alpha: slope in the negative part of the relu. Defaults to 0.01.
#' @returns the compiled neural network model
build_model <- function(
    input_dim=5,
    hidden_units=c(20,10),
    initial_lr = 0.01,
    use_tanh=FALSE,
    leaky_relu_alpha=0.1){
  model <- -1 # reset the model variable.
  # definition of the model
  if(use_tanh){
    model <- keras_model_sequential() %>%
      layer_dense(units = hidden_units[1], input_shape = input_dim, activation='tanh') %>%
      layer_dense(units = hidden_units[2], activation='tanh') %>%
      layer_dense(units = 1, activation = 'sigmoid') # sigmoid activation, for binary classification
  } else {
    model <- keras_model_sequential() %>%
      layer_dense(units = hidden_units[1], input_shape = input_dim) %>%
      layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
      layer_dense(units = hidden_units[2]) %>%
      layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
      layer_dense(units = 1, activation = 'sigmoid') # sigmoid activation, for binary classification
  }
  # compile the model, define optimizer, loss and metrics
  model %>% compile(
    optimizer = keras$optimizers$Adam(learning_rate=initial_lr),
    loss = keras$losses$BinaryCrossentropy(),
    metrics = keras$metrics$BinaryAccuracy()
  )
  return(model)
}

#' function to train the model
#' @param model: The NN model to train
#' @param x_train: training data
#' @param y_train: labels for the training data
#' @param lr_schedule: Function for learning rate scheduling. Defaults to lr_schedule_fun.
#' @param num_epochs: Integer, defaults to 500.
#' Number of epochs for which the model should be trained.
#' @param batch_size: Integer, defaults to 128.
#' Batch size to use to compute gradients during training
#' @param val_split: Number between 0 and 1. Fraction of data used for validation.
#' @param verbose: 0 or 1, defaults to 1.
#' If 1, an output is printed after every epoch of the training process.
#' If 0, no output is printed.
#' @returns model: the trained model
train_model <- function(
    model,
    x_train,
    y_train,
    lr_schedule = lr_schedule_fun,
    num_epochs = 300,
    batch_size=128,
    val_split = 0.2,
    verbose=1
){
  # create the keras object required for learning rate scheduling
  lr_scheduler <- callback_learning_rate_scheduler(schedule = lr_schedule)

  # train the model
  history <- model %>% fit(
    x_train, y_train,
    epochs = num_epochs,
    batch_size = batch_size,
    validation_split=val_split, # use 20 percent of training data as validation data
    verbose = 1, # 0 for slightly faster training (no output), 1 to observe progress while training
    callbacks=list(lr_scheduler)
  )
  return(model)
}

#' Compute the g values as defined in the thesis
#' @param model: Neural Network for binary classification
#' @param orig_data: The observed data.
#' @param nu: Defaults to 1. Fraction of noise samples to real samples that
#' was used during training of the neural network.
#' @returns A vector containing the values g_t (see thesis)
compute_gvals <- function(
    model,
    orig_data,
    nu=1){
  predictions <- model %>% predict(orig_data)
  length(predictions[predictions> 0.5]) / length(predictions)
  # here a classifier p(u) was trained on the data. to get the values g(u),
  # as defined in the thesis, g(u) = log(\nu \cdot p(u)/(1-p(u))) needs to be computed,
  # where nu = T_n/T_c, the fraction of number of noise samples to observed samples
  g_vals <- log(nu * predictions / (1-predictions)) # G(u, \eta) from the thesis
  return(g_vals)
}

#' performs a quantile regression
#' @param g_vals: Vector. Values g_t from the thesis
#' @param orig_data: The observed data.
#' @param family_set_name: Defaults to "onepar".
#' Argument family_set in the D-vine Quantile regression function.
#' Possible choices include "onepar" and "parametric".
#' @param bottom_quantile_levels: Defaults to c(0.05,0.1). Vector of lower quantiles,
#' for which it is checked, whether the conditioned quantile estimates are >0.
#' Tests whether the alternative model is better
#' @param top_quantile_levels: Defaults to c(0.05,0.1). Vector of upper quantiles,
#' for which it is checked, whether the conditioned quantile estimates are <0.
#' Tests whether the simplified model is better
#' @returns: A List of 2 vectors, the first containing the number of samples, for which
#' q(bottom_quantile_levels) > 0 holds, the second containing the number of samples, for
#' which q(top_quantile_levels) < 0 holds.
perform_quant_reg <- function(
    g_vals,
    orig_data,
    family_set_name = "onepar",
    bottom_quantile_levels = c(0.05,0.1),
    top_quantile_levels = c(0.9,0.95)){
  orig_data <- data.frame(orig_data)
  qreg_data <- cbind(g_vals, orig_data)
  q_reg <- vinereg::vinereg(g_vals ~ . ,family_set=family_set_name, data=qreg_data)
  bottom_quantiles <- predict(q_reg, newdata=orig_data, alpha=bottom_quantile_levels)
  top_quantiles <- predict(q_reg, newdata=orig_data, alpha=top_quantile_levels)
  alternative_better <- rep(0,length(bottom_quantile_levels))
  simp_better <- rep(0,length(top_quantile_levels))
  for(i in 1:nrow(orig_data)){
    for(j in 1:length(bottom_quantile_levels)){
      if(bottom_quantiles[i,j] > 0){
        alternative_better[j] <- alternative_better[j] + 1
      }
    }
    for(j in 1:length(top_quantile_levels)){
        if(top_quantiles[i,j] < 0){
          simp_better[j] <- simp_better[j] + 1
        }
    }
  }
  output <- list(alternative_better,
                 simp_better)
  return(output)

}


# Example Usage:
# split_output <- train_test_split(orig_data, simplified_data)
# x_train <- split_output[[1]]
# x_test <- split_output[[2]]
# y_train <- split_output[[3]]
# y_test <- split_output[[4]]
# model <- build_model(input_dim=ncol(x_train))
# model %>% summary
# model <- train_model(model, x_train, y_train, num_epochs=10)
# model %>% evaluate(x_test, y_test)
# g_vals <- compute_gvals(model, orig_data)
# output <- perform_quant_reg(g_vals, orig_data)
# print(output)
