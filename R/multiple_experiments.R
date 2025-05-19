# Script to run many simulations
source("R/simulate_non_simplified_vine.R")
source("R/classifier_methods.R")

num_samples <- list(1000,10000) #2 different number of samples
nus <- list(10,4) # n_{noise} / n_{true}. Higher, if only 1000 samples given.
dims <- list(3,4,5) # dimensions
# each list element contains first the lower, then the upper limit
tau_limits <- list(c(0.001,0.6), c(0.001,0.9))

# Define 3d parameters
struct_mat_3d <- matrix(c(1,1,1,
                         2,2,0,
                         3,0,0), ncol=3, byrow=TRUE)
families_3d <- list(list("frank", "gaussian"), list("frank"))
params_3d <- list(c(ktau_to_par(family=families_3d[[1]][[2]], tau=-0.2)),
                  c(ktau_to_par(family=families_3d[[1]][[2]], tau=0.3)))
rotations_3d <- list(list(0,0),0)
# Define 4d parameters
struct_mat_4d <- matrix(c(3,3,3,3,
                         2,2,2,0,
                         4,4,0,0,
                         1,0,0,0), ncol=4, byrow=TRUE)
families_4d <- list(list("frank", "gaussian","gaussian"),
                   list("gaussian","frank"),
                   list("frank"))
params_4d <- list(c(ktau_to_par(family=families_4d[[1]][[2]], tau=-0.2)),
                 c(ktau_to_par(family=families_4d[[1]][[2]], tau=0.3)),
                 c(ktau_to_par(family=families_4d[[1]][[3]], tau=0.7)))
rotations_4d <- list(list(0,0,0),list(0,0), list(0))
# Define 5d parameters
struct_mat_5d <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
families_5d <- list(list("frank", "gaussian","gaussian","frank"),
                    list("frank","gaussian","gaussian"),
                    list("gaussian", "frank"),
                    list("gaussian"))
params_5d <- list(c(ktau_to_par(family=families_5d[[1]][[1]], tau=-0.2)),
                  c(ktau_to_par(family=families_5d[[1]][[2]], tau=0.3)),
                  c(ktau_to_par(family=families_5d[[1]][[3]], tau=-0.1)),
                  c(ktau_to_par(family=families_5d[[1]][[4]], tau=0.1)))
rotations_5d <- list(list(0,0,0,0),list(0,0,0), list(0,0), list(0))
# Summarize parameters by dimension
struct_mats <- list(struct_mat_3d, struct_mat_4d, struct_mat_5d)
families <- list(families_3d, families_4d, families_5d)
initial_params <- list(params_3d, params_4d, params_5d)
rotations <- list(rotations_3d, rotations_4d, rotations_5d)
# Quantile Levels to calculate:
lower_q_levs <- seq(0.01,0.2,0.01)
upper_q_levs <- seq(0.8,0.99,0.01)
# Run big simulation
total_runs <- 0
# Initialize an empty list to store the results of each inner loop iteration
all_results <- list()
for(sample_idx in 1:length(num_samples)){
  # Check, if any of the dimensions are not implemented. If so stop execution and print.
  incorrect_dimensions <- setdiff(unlist(dims), c(3,4,5))
  if (length(incorrect_dimensions) > 0) {
    print("The following dimensions are not implemented: ")
    print(incorrect_dimensions)
    print("Please remove these dimensions, then try to execute this function again.")
    break
  }
  nu_var <- nus[[sample_idx]]
  for(tau_idx in 1:length(tau_limits)){
    tau_lower=tau_limits[[tau_idx]][1]
    tau_upper=tau_limits[[tau_idx]][2]
    for(dim_idx in 1:length(dims)){
      # initialize param_cond_funcs (overwritten below)
      param_cond_func_var <- -1
      # Define the conditional parameter functions
      if(dims[[dim_idx]] ==3){
        param_cond_func_3d_lin <- list(
          list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_func_3d_non_lin <- list(
          list(u_to_param_non_lin(c(1,1),
                                  tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_func_3d_non_lin_cubic <- list(
          list(u_to_param_non_lin_cubic(c(1,1.2,0.5),
                                        tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_func_var <- list(
          param_cond_func_3d_lin,
          param_cond_func_3d_non_lin,
          param_cond_func_3d_non_lin_cubic)
      } else if(dims[[dim_idx]] ==4){
        param_cond_funcs_4d_lin <- list(
          list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_funcs_4d_non_lin <- list(
          list(
            u_to_param_non_lin(c(1,1.4),
                               tau_lower=tau_lower, tau_upper=tau_upper),
            u_to_param_non_lin(c(0.6,0.8),
                               tau_lower=tau_lower, tau_upper=tau_upper)),
          list(
            u_to_param_non_lin(c(0.4,0.8, -0.3,0.7,1.3),
                               tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_funcs_4d_non_lin_cubic <- list(
          list(u_to_param_non_lin_cubic(c(1,0.5,0.7),
                                        tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_non_lin_cubic(c(0.4,0.4,0.3),
                                        tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_non_lin_cubic(c(0.4,0.7, 0.2,-0.4,0.8,-0.9,0.6,0.2,0.4),
                                        tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_func_var <- list(
          param_cond_funcs_4d_lin,
          param_cond_funcs_4d_non_lin,
          param_cond_funcs_4d_non_lin_cubic)
      } else if(dims[[dim_idx]] ==5){
        param_cond_funcs_5d_lin <- list(
          list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_linear(c(0.2,0.5,0.3), tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_funcs_5d_non_lin <- list(
          list(
            u_to_param_non_lin(c(1,1),
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
        param_cond_funcs_5d_non_lin_cubic <- list(
          list(u_to_param_non_lin_cubic(c(1,1,0.5),
                                        tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_non_lin_cubic(c(1,0.7,0.3),
                                        tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_non_lin_cubic(c(1,2,0.1),
                                        tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_non_lin_cubic(c(0.7,0.3, 0.9,-0.4,0.8,0.1,0.2,0.3,0.4),
                                        tau_lower=tau_lower, tau_upper=tau_upper),
               u_to_param_non_lin_cubic(c(0.4,0.6,1,0.7,1.8,0.01,0.04,0.07,0.03),
                                        tau_lower=tau_lower, tau_upper=tau_upper)),
          list(u_to_param_non_lin_cubic(c(0.2, 0.5, 0.3,
                                          1,   1.4, 1  ,0.4,0.7,0.8,
                                          0.01,0.04,0.07,0.03,0.02,0.05,0.07,0.04,0.02,0.1),
                                        tau_lower=tau_lower, tau_upper=tau_upper)))
        param_cond_func_var <- list(
          param_cond_funcs_5d_lin,
          param_cond_funcs_5d_non_lin,
          param_cond_funcs_5d_non_lin_cubic)
      }
      # Define the other necessary parameters for simulation
      struct_mat_var <- struct_mats[[dim_idx]]
      families_var <- families[[dim_idx]]
      params_var <- initial_params[[dim_idx]]
      rotations_var <- rotations[[dim_idx]]
      for(par_idx in 1:length(param_cond_func_var)){
        total_runs <- total_runs+1
        print(paste0("Executing... ",
                     total_runs,
                     " of a total of ",
                     print(length(dims)*length(num_samples)*length(tau_limits)*length(param_cond_func_var))))
        non_simp_data <- simulate_non_simp_parallel(n_samples = num_samples[[sample_idx]],
                                             struct = struct_mat_var,
                                             families=families_var,
                                             params = params_var,
                                             param_cond_funcs = param_cond_func_var[[par_idx]],
                                             rotations = rotations_var)
        orig_data <- as.matrix(non_simp_data)
        orig_data <- unname(orig_data)
        fitted_vine <-vinecop(orig_data,family_set="onepar")
        # number of samples to create from the fitted vine
        num_samples <- nu_var * nrow(orig_data)
        # simulate from the simplified vine, to train a classifier
        simplified_samples <- rvinecop(num_samples, fitted_vine)
        split_output <- train_test_split(orig_data, simplified_samples, train_perc=0.8, stratified=TRUE)
        x_train <- split_output[[1]]
        x_test <- split_output[[2]]
        y_train <- split_output[[3]]
        y_test <- split_output[[4]]
        #sum(y_train) /length(y_train)
        # delete model from previous iteration
        model <- 0
        # define new model
        model <- build_model(input_dim=ncol(x_train), use_tanh=FALSE)
        # initial evaluation (to test implementation)
        pre_training_eval <- model %>% evaluate(x_train, y_train)
        model <- train_model(model, x_train, y_train, num_epochs=100)
        test_set_eval <- model %>% evaluate(x_test, y_test)
        g_vals <- compute_gvals(model, orig_data, nu=nu_var)
        output <- perform_quant_reg(g_vals,
                                    orig_data,
                                    bottom_quantile_levels=lower_q_levs,
                                    top_quantile_levels = upper_q_levs)
        # Store the results of the current iteration
        current_result <- c(num_samples[[sample_idx]],
                            dims[[dim_idx]],
                            tau_limits[[tau_idx]],
                            par_idx,
                            output[[1]],
                            output[[2]])
        # append to the results list, as a list, to keep the rows separated
        all_results <- append(all_results, list(current_result))

      }
    }
  }
}
print(total_runs)
print(length(dims)*length(num_samples)*length(tau_limits)*3)
# Summarize the results in a dataframe
results_df <- as.data.frame(do.call(rbind, all_results))
# Give better names to the columns of the dataframe
lower_quantile_names <- paste0("q(", lower_q_levs, ")>0")
upper_quantile_names <- paste0("q(", upper_q_levs, ")<0")
colnames(results_df) <- c("num_samples",
                          "dim",
                          "tau_lower",
                          "tau_upper",
                          "param_cond_func_idx",
                          lower_quantile_names,
                          upper_quantile_names)
# print the result dataframe.
print(head(results_df))
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
# Construct the file name with the date
results_filename <- paste0("results/", current_date, ".csv")
write.csv(results_df, file = results_filename, row.names = FALSE)
print(paste("Data saved to:", results_filename))


orig_data <- as.matrix(u_test)
orig_data <- unname(orig_data)
nu = 10
fitted_vine <-vinecop(orig_data,family_set="onepar")
# simulate from the simplified vine, to train a classifier later.
num_samples <- nu * nrow(orig_data) # number of samples to create from the fitted vine
simplified_samples <- rvinecop(num_samples, fitted_vine)
split_output <- train_test_split(orig_data, simplified_samples, train_perc=0.8, stratified=TRUE)
x_train <- split_output[[1]]
x_test <- split_output[[2]]
y_train <- split_output[[3]]
y_test <- split_output[[4]]
sum(y_train) /length(y_train)
model <- 0
model <- build_model(input_dim=ncol(x_train), use_tanh=TRUE)
model %>% summary
# initial evaluation (to test implementation)
pre_training_eval <- model %>% evaluate(x_train, y_train)
pre_training_eval
model <- train_model(model, x_train, y_train, num_epochs=100)
test_set_eval <- model %>% evaluate(x_test, y_test)
test_set_eval[["loss"]]
test_set_eval[["binary_accuracy"]]
g_vals <- compute_gvals(model, orig_data, nu=nu)
lower_q_levs <- seq(0.01,0.2,0.01)
lower_q_levs
upper_q_levs <- seq(0.8,0.99,0.01)
upper_q_levs
output <- perform_quant_reg(g_vals,
                            orig_data,
                            bottom_quantile_levels=lower_q_levs,
                            top_quantile_levels = upper_q_levs)
print(output)
print(output[[1]] / nrow(orig_data))
plot(lower_q_levs, output[[1]])
plot(upper_q_levs, output[[2]])
