# Script to run many simulations
source("R/simulate_non_simplified_vine.R")
source("R/classifier_methods.R")

num_samples <- c(1000,10000)
dims <- c(3,4,5)
tau_vals_lower <- c(0.001,0.001)
tau_vals_upper <- c(0.5,0.9)
# TODO: write script that loops over all the parameters and automatically performs the training.

struct_mat <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
#plot(rvine_matrix(struct_mat),1:4)
family_test <- list(list("frank", "clayton","gaussian","frank"),
                    list("frank","gaussian","joe"),
                    list("gaussian", "gumbel"),
                    list("gaussian"))
params_test <- list(c(ktau_to_par(family=family_test[[1]][[1]], tau=-0.2)),
                    c(ktau_to_par(family=family_test[[1]][[2]], tau=0.3)),
                    c(ktau_to_par(family=family_test[[1]][[3]], tau=-0.1)),
                    c(ktau_to_par(family=family_test[[1]][[4]], tau=0.1)))
tau_lower = 0.001
tau_upper = 0.92
param_cond_funcs_test <- list(list(u_to_param_non_lin(c(1,1),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin(c(1,0.7),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin(c(1,2),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_non_lin(c(0.7,0.3, 0.9,-0.4,0.8),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin(c(0.4,0.6,1,0.7,2),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_non_lin(c(0.2,0.5,0.3,1,2,1,0.4,0.7,0.8),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)))
u_test <- simulate_non_simp_parallel(n_samples = 5000,
                                     struct = struct_mat,
                                     families=family_test,
                                     params = params_test,
                                     param_cond_funcs = param_cond_funcs_test,
                                     rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
#pairs_copula_data(u_test)

orig_data <- as.matrix(u_test)
orig_data <- unname(orig_data)
nu = 1
fitted_vine <-vinecop(orig_data,family_set="onepar")
# simulate from the simplified vine, to train a classifier later.
num_samples <- nu * nrow(orig_data) # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
split_output <- train_test_split(orig_data, simplified_samples)
x_train <- split_output[[1]]
x_test <- split_output[[2]]
y_train <- split_output[[3]]
y_test <- split_output[[4]]

model <- build_model(input_dim=ncol(x_train))
model %>% summary
model <- train_model(model, x_train, y_train, num_epochs=100)
model %>% evaluate(x_test, y_test)
g_vals <- compute_gvals(model, orig_data)
output <- perform_quant_reg(g_vals, orig_data)
print(output)
print(output[[1]] / nrow(orig_data))
