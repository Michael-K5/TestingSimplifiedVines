# Fit a simplified vine copula and train a classifier on the data
library(rvinecopulib)
library(keras)
library(tensorflow)
# load the original data
csv_filename <- "data/non_simplified_sim_2025-04-15.csv"
orig_data <- read.csv(csv_filename)
num_rows <- nrow(orig_data)
# cast labels as matrix (otherwise this does not work with tensorflow)
labels <- as.matrix(rep(0L, num_rows))
#pairs_copula_data(orig_data)

# fit a simplified vine
fitted_vine <-vinecop(orig_data,family_set="onepar")
print.data.frame(summary(fitted_vine),digit=2)
num_samples <- num_rows # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
#pairs_copula_data(simplified_samples)
simplified_labels <- as.matrix(rep(1L, num_samples))

classifier_data <- rbind(orig_data, simplified_samples)
classifier_labels <- rbind(labels, simplified_labels)

# train test split
train_perc <- 0.8
num_train <- floor(train_perc * nrow(classifier_data))
train_idx <- sample(nrow(classifier_data),num_train)
x_train <- classifier_data[train_idx,]
y_train <- as.matrix(classifier_labels[train_idx,])
#y_train <- tf$keras$utils$to_categorical(classifier_labels[train_idx], num_classes = 2L)
x_test <- classifier_data[-train_idx,]
#y_test <- tf$keras$utils$to_categorical(classifier_labels[-train_idx], num_classes = 2L)
y_test <- as.matrix(classifier_labels[-train_idx,])
#length(y_test)

Sequential <- tf$keras$models$Sequential
Dense <- tf$keras$layers$Dense
Input <- tf$keras$layers$Input
Dropout <- tf$keras$layers$Dropout

# Initialize the model
model <- Sequential()
# Add layers to the model
model$add(Dense(units = 20L, activation = 'tanh', input_shape=list(as.integer(ncol(x_train)))))
#model$add(Dropout(rate=0.3))
model$add(Dense(units = 10L, activation = 'tanh'))
#model$add(Dropout(rate=0.3))
model$add(Dense(units = 1L, activation = 'sigmoid'))

# compile the model, define loss, optimizer and metrics
model$compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam',
  metrics = list('accuracy')
)

model$summary()

# fit the model to the data
history <- model$fit(
  x_train, y_train,
  epochs = 500L, batch_size = 100L,
  validation_split = 0.2
)

model$evaluate(x_test,y_test)

plot(history$history$val_accuracy, ylab="validation accuracy", xlab="epoch")
