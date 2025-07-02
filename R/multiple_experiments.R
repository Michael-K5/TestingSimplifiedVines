# Script to run many simulations
source("R/simulate_non_simplified_vine.R")
source("R/classifier_methods.R")
# Parameters for the simulation.
# 2 different number of samples.
num_samples <- list(1000,10000)
# fraction of number of noise to true samples. Higher, if only 1000 samples given.
nus <- list(10,4)
# dimensions of the simulated data.
dims <- list(3,4,5)
# each list element contains first the lower, then the upper limit
tau_limits <- list(c(0.001,0.3), c(0.001,0.6), c(0.001,0.9))

# Define 3d parameters
struct_mat_3d <- matrix(c(1,1,1,
                         2,2,0,
                         3,0,0), ncol=3, byrow=TRUE)
families_3d <- list(list("frank", "gaussian"), list("frank"))
params_3d <- list(c(ktau_to_par(family=families_3d[[1]][[1]], tau=0.2)),
                  c(ktau_to_par(family=families_3d[[1]][[2]], tau=0.4)))
rotations_3d <- list(list(0,0),0)
# Define 4d parameters
struct_mat_4d <- matrix(c(3,3,3,3,
                         2,2,2,0,
                         4,4,0,0,
                         1,0,0,0), ncol=4, byrow=TRUE)
families_4d <- list(list("frank", "gaussian","gaussian"),
                   list("gaussian","frank"),
                   list("frank"))
params_4d <- list(c(ktau_to_par(family=families_4d[[1]][[1]], tau=0.2)),
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
params_5d <- list(c(ktau_to_par(family=families_5d[[1]][[1]], tau=0.2)),
                  c(ktau_to_par(family=families_5d[[1]][[2]], tau=0.3)),
                  c(ktau_to_par(family=families_5d[[1]][[3]], tau=0.4)),
                  c(ktau_to_par(family=families_5d[[1]][[4]], tau=0.1)))
rotations_5d <- list(list(0,0,0,0),list(0,0,0), list(0,0), list(0))
# Summarize parameters by dimension
# struct_mats <- list(struct_mat_3d, struct_mat_4d, struct_mat_5d)
# families <- list(families_3d, families_4d, families_5d)
# initial_params <- list(params_3d, params_4d, params_5d)
# rotations <- list(rotations_3d, rotations_4d, rotations_5d)

#struct_mats <- list(struct_mat_4d)
#families <- list(families_4d)
#initial_params <- list(params_4d)
#rotations <- list(rotations_4d)

#' Run the simulations
#' @param num_samples: list of integers, that determine for what different
#' number of samples drawn from a non-simplified vine copula the experiments should be run.
#' @param nus: list of numbers: Determines how many noise (simplified) samples should be simulated.
#' The i-th entry contains the values nu = n_{noise}/n_{true} for
#' the number of samples specified in num_samples[[i]]
#' @param dims: list of integers, that determine for what dimensions to run the tests
#' (currently implemented for 3, 4 and 5)
#' @param tau_limits: list of 2 dimensional vectors: For each of those vectors,
#' the first entry determines, what the lower threshold of kendalls tau is,
#' the second entry determines what the upper threshold of kendalls tau is in the
#' conditional copulas
#' @param struct_mats: List of regular vine matrices: Determine what structure to use.
#' The entry in position i needs to have the same dimension as dims[[i]].
#' @param families: List of (list of list of string): The i-th element
#' contains the copula families corresponding to the copulas defined with struct_mats[[i]]
#' @param initial_params: List of vectors: The i-th entry contains the parameters
#' for the i-th unconditioned copula in the first tree, defined by struct_mats[[i]]
#' @param rotations: List of (list of list of int): The i-th element contains the
#' rotations of the copulas specified in struct_mats[[i]]
#' @param lower_q_levs: Vector of float: defaults to seq(0.01,0.25,0.01). The lower
#' quantiles that should be estimated in the quantile regression
#' @param upper_q_levs: Vector of float: defaults to seq(0.75,0.99,0.01). The upper
#' quantiles that should be estimated in the quantile regression
#' @param int_bottom_threshold: Float, defaults to 0.6: If the monte carlo integral of the
#' non-parametric model is less than int_bottom_threshold, then the same combination of parameters
#' is run again (at most max_retries times), to try to find a more appropriate model.
#' @param int_top_threshold: Float, defaults to 1.8: If the monte carlo integral of the
#' non-parametric model is greater than int_top_threshold, then the same combination of parameters
#' is run again (at most max_retries times), to try to find a more appropriate model.
#' @param max_retries: Integer, defaults to 3: How often to retry evaluating the model if
#' the monte carlo integral value is not close enough to 1 (as measured by int_dev_threshold)
#' @param filename: string, defaults to "": If not "", then the results are written to
#' a csv file, with name specified as filename
#' @returns result_df: A Dataframe with the results of the experiments.
run_simulations <- function(num_samples,
                            dims,
                            nus,
                            tau_limits,
                            struct_mats,
                            families,
                            initial_params,
                            rotations,
                            lower_q_levs = seq(0.01,0.25,0.01),
                            upper_q_levs = seq(0.75,0.99,0.01),
                            bottom_q_NN_train = c(0.01,0.05,0.1,0.15,0.2,0.25),
                            top_q_NN_train=c(0.75,0.8,0.85,0.9,0.95,0.99),
                            int_bottom_threshold = 0.6,
                            int_top_threshold = 1.9,
                            max_retries = 3,
                            filename=""){
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
            list(u_to_param_quadratic(c(1,1),
                                    tau_lower=tau_lower, tau_upper=tau_upper)))
          param_cond_func_3d_non_lin_cubic <- list(
            list(u_to_param_cubic(c(1,1.2,0.5),
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
              u_to_param_quadratic(c(1,1.4),
                                 tau_lower=tau_lower, tau_upper=tau_upper),
              u_to_param_quadratic(c(1.2,0.9),
                                 tau_lower=tau_lower, tau_upper=tau_upper)),
            list(
              u_to_param_quadratic(c(-0.7,1.5,0.9,-1.2,1.5),
                                 tau_lower=tau_lower, tau_upper=tau_upper)))
          param_cond_funcs_4d_non_lin_cubic <- list(
            list(u_to_param_cubic(c(0.7,1.2,-0.8),
                                          tau_lower=tau_lower, tau_upper=tau_upper),
                 u_to_param_cubic(c(1.4,-1.8,1.0),
                                          tau_lower=tau_lower, tau_upper=tau_upper)),
            list(u_to_param_cubic(c(0.7,1.8, -0.9,1.4,0.8,1.1,-1.2,-1.3,1.7),
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
              u_to_param_quadratic(c(1,1),
                                 tau_lower=tau_lower, tau_upper=tau_upper),
              u_to_param_quadratic(c(1,0.7),
                                 tau_lower=tau_lower, tau_upper=tau_upper),
              u_to_param_quadratic(c(1,2),
                                 tau_lower=tau_lower, tau_upper=tau_upper)),
            list(u_to_param_quadratic(c(0.7,0.5, 0.9,-1.2,1.5),
                                    tau_lower=tau_lower, tau_upper=tau_upper),
                 u_to_param_quadratic(c(0.4,-0.8,1.0,0.7,2.0),
                                    tau_lower=tau_lower, tau_upper=tau_upper)),
            list(u_to_param_quadratic(c(0.7,0.4,-0.9,1.3,0.75,1.3,-0.6,0.5,-1.1),
                                    tau_lower=tau_lower, tau_upper=tau_upper)))
          param_cond_funcs_5d_non_lin_cubic <- list(
            list(u_to_param_cubic(c(1.0,1.0,0.5),
                                          tau_lower=tau_lower, tau_upper=tau_upper),
                 u_to_param_cubic(c(1.0,0.7,0.3),
                                          tau_lower=tau_lower, tau_upper=tau_upper),
                 u_to_param_cubic(c(1.0,2.0,1.1),
                                          tau_lower=tau_lower, tau_upper=tau_upper)),
            list(u_to_param_cubic(c(0.7,0.8, -0.9,0.4,0.8,1.1,-1.2,-1.3,1.0),
                                          tau_lower=tau_lower, tau_upper=tau_upper),
                 u_to_param_cubic(c(0.4,0.6,1,0.7,2, -1.1,1.2,-0.9,-0.3),
                                          tau_lower=tau_lower, tau_upper=tau_upper)),
            list(u_to_param_cubic(c(0.7, 0.5, -1.3,
                                    1,   1.4, 1  ,-1.4,0.7,-0.8,
                                    1.1,-1.4,0.7,-0.3,0.4,
                                    0.8, -0.7,1.4,-1.2,-0.6),
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
                       length(dims)*length(num_samples)*length(tau_limits)*length(param_cond_func_var)))
          non_simp_data <- simulate_non_simplified(n_samples = num_samples[[sample_idx]],
                                                      struct = struct_mat_var,
                                                      families=families_var,
                                                      params = params_var,
                                                      param_cond_funcs = param_cond_func_var[[par_idx]],
                                                      rotations = rotations_var)
          orig_data <- as.matrix(non_simp_data)
          orig_data <- unname(orig_data)
          fitted_vine <-vinecop(orig_data,family_set="parametric")
          # number of samples to create from the fitted vine
          n_noise <- nu_var * nrow(orig_data)
          # simulate from the simplified vine, to train a classifier
          simplified_samples <- rvinecop(n_noise, fitted_vine)
          split_output <- train_test_split(orig_data=orig_data,
                                           simplified_data=simplified_samples,
                                           train_perc=0.8,
                                           stratified=TRUE)
          x_train <- split_output[[1]]
          x_test <- split_output[[2]]
          y_train <- split_output[[3]]
          y_test <- split_output[[4]]
          int_val <- int_bottom_threshold - 1
          retry_num <- 0
          model <- 0
          train_set_eval <- -1
          test_set_eval <- -1
          # Loop to discard a certain value, if the integral deviates too much from 1
          while((int_val < int_bottom_threshold || int_val > int_top_threshold)
                && retry_num < max_retries){
            retry_num <- retry_num +1
            # define new model
            model <- build_model(input_dim=ncol(x_train), use_tanh=FALSE)
            train_model_output <- train_model(model, x_train, y_train, num_epochs=200)
            model <- train_model_output[[1]]
            test_set_eval <- model %>% evaluate(x_test, y_test)
            train_set_eval <- model %>% evaluate(x_train, y_train)
            int_val <- compute_integral(model=model,
                                        fitted_vine = fitted_vine,
                                        nu=nu_var,
                                        data_dim_if_unif = ncol(orig_data),
                                        n_samples=20000,
                                        user_info=FALSE)
          }
          log_lik_true_train <- log_likelihood_non_simplified(
            u_data=x_train,
            struct = struct_mat_var,
            families=families_var,
            params = params_var,
            param_cond_funcs = param_cond_func_var[[par_idx]],
            rotations = rotations_var,
            return_vector=FALSE)
          log_lik_true_test <- log_likelihood_non_simplified(
            u_data=x_test,
            struct = struct_mat_var,
            families=families_var,
            params = params_var,
            param_cond_funcs = param_cond_func_var[[par_idx]],
            rotations = rotations_var,
            return_vector=FALSE)
          # evaluate log likelihoods of original data (on train and test set fractions of orig_data)
          log_lik_simp_train_orig <- sum(log(dvinecop(x_train[y_train==1,], fitted_vine)))
          log_lik_NN_train_orig <- sum(log(non_param_cop(
            model=model,
            fitted_vine=fitted_vine,
            obs=x_train[y_train==1,],
            nu=nu_var)))
          log_lik_simp_test_orig <- sum(log(dvinecop(x_test[y_test==1,], fitted_vine)))
          log_lik_NN_test_orig <- sum(log(non_param_cop(
            model=model,
            fitted_vine=fitted_vine,
            obs=x_test[y_test==1,],
            nu=nu_var)))
          # get number of model parameters and degrees of freedom for the log likelihood ratio test
          NN_num_params <- count_NN_params(weights=model$weights)
          simp_cop_num_params <- get_num_cop_params(fitted_vine)
          deg_free <- NN_num_params # 20250702 changed! removed -simp_cop_num_params
          # log likelihood ratio test statistic
          LRT_stat_train <- -2*(as.numeric(log_lik_simp_train_orig) - as.numeric(log_lik_NN_train_orig))
          LRT_stat_test <- -2*(as.numeric(log_lik_simp_test_orig) - as.numeric(log_lik_NN_test_orig))
          LRT_stat_overall <- LRT_stat_train + LRT_stat_test
          # compute p(chisq(df=deg_free) > LRT_stat) (for train and test set)
          p_value_train <- pchisq(LRT_stat_train, df = deg_free, lower.tail = FALSE)
          p_value_test <- pchisq(LRT_stat_test, df = deg_free, lower.tail = FALSE)
          p_value_overall <- pchisq(LRT_stat_overall, df=deg_free, lower.tail=FALSE)
          # compute AIC and BIC for the NN and simplified model
          AIC_NN <- 2*(NN_num_params + simp_cop_num_params) - 2* (log_lik_NN_train_orig + log_lik_NN_test_orig)
          AIC_simp <- 2*simp_cop_num_params - 2*(log_lik_simp_train_orig + log_lik_simp_test_orig)
          BIC_NN <- log(nrow(orig_data))*(NN_num_params + simp_cop_num_params) - 2* (log_lik_NN_train_orig + log_lik_NN_test_orig)
          BIC_simp <- log(nrow(orig_data))*simp_cop_num_params - 2*(log_lik_simp_train_orig + log_lik_simp_test_orig)
          # End likelihood ratio test
          # paired wilcoxon and t test
          log_cor_facs <- log(correction_factors(model, obs=orig_data, nu=nu_var))
          wilcox_results <- wilcox.test(
            log_cor_facs,
            alternative="two.sided",
            mu=0,
            exact=FALSE,
            conf.int=TRUE,
            conf.level=0.95
          )
          wilcox_test_statistic <- wilcox_results$statistic
          wilcox_p_value <- wilcox_results$p.value
          wilcox_conf_int_bottom <- wilcox_results$conf.int[1]
          wilcox_conf_int_top <- wilcox_results$conf.int[2]
          t_test_results <- t.test(
            log_cor_facs,
            alternative="two.sided",
            mu=0
          )
          t_test_statistic <- t_test_results$statistic
          t_test_p_value <- t_test_results$p.value
          t_test_conf_int_bottom <- t_test_results$conf.int[1]
          t_test_conf_int_top <- t_test_results$conf.int[2]
          # End Wilcoxon and t test
          # Perform D-Vine Quantile regression (parametric)
          g_vals <- compute_gvals(model, orig_data, nu=nu_var)
          qreg_par_output <- perform_quant_reg(orig_data,
                                      g_vals,
                                      family_set_name="onepar",
                                      bottom_quantile_levels=lower_q_levs,
                                      top_quantile_levels = upper_q_levs)
          # Perform D-Vine Quantile regression (non-parametric)
          qreg_nonpar_output <- perform_quant_reg(orig_data,
                                      g_vals,
                                      family_set_name="nonparametric",
                                      bottom_quantile_levels=lower_q_levs,
                                      top_quantile_levels = upper_q_levs)
          # Linear Quantile Regression
          lin_qreg_output <- perform_linear_quant_reg(
            orig_data,
            g_vals,
            bottom_quantiles_lin =lower_q_levs,
            top_quantiles_lin =upper_q_levs,
            train_perc=0.9)
          # MCQRNN
          mcqrnn_output <- perform_quant_reg_mcqrnn(
            orig_data=orig_data,
            g_vals=g_vals,
            bottom_q_NN_train = bottom_q_NN_train,
            bottom_q_NN_predict=lower_q_levs,
            top_q_NN_train=top_q_NN_train,
            top_q_NN_predict = upper_q_levs,
            num_hidden = 10,
            train_perc=0.9,
            user_info=FALSE
          )
          # Store the results of the current iteration
          current_result <- c(num_samples[[sample_idx]],
                              nus[[sample_idx]],
                              dims[[dim_idx]],
                              tau_limits[[tau_idx]], # this is a vector with 2 elements
                              par_idx,
                              qreg_par_output[[1]],
                              qreg_par_output[[2]],
                              qreg_par_output[[3]],
                              qreg_nonpar_output[[1]],
                              qreg_nonpar_output[[2]],
                              qreg_nonpar_output[[3]],
                              lin_qreg_output[[1]],
                              lin_qreg_output[[2]],
                              lin_qreg_output[[3]],
                              lin_qreg_output[[4]],
                              mcqrnn_output[[1]],
                              mcqrnn_output[[2]],
                              mcqrnn_output[[3]],
                              mcqrnn_output[[4]],
                              train_set_eval[[1]],
                              train_set_eval[[2]],
                              test_set_eval[[1]],
                              test_set_eval[[2]],
                              int_val,
                              log_lik_simp_train_orig,
                              log_lik_NN_train_orig,
                              log_lik_true_train,
                              log_lik_simp_test_orig,
                              log_lik_NN_test_orig,
                              log_lik_true_test,
                              NN_num_params,
                              simp_cop_num_params,
                              deg_free,
                              LRT_stat_train,
                              LRT_stat_test,
                              LRT_stat_overall,
                              p_value_train,
                              p_value_test,
                              p_value_overall,
                              AIC_NN,
                              AIC_simp,
                              BIC_NN,
                              BIC_simp,
                              wilcox_test_statistic,
                              wilcox_p_value,
                              wilcox_conf_int_bottom,
                              wilcox_conf_int_top,
                              t_test_statistic,
                              t_test_p_value,
                              t_test_conf_int_bottom,
                              t_test_conf_int_top
                              )
          # append to the results list, as a list, to keep the rows separated
          all_results <- append(all_results, list(current_result))
        }
      }
    }
  }
  # Summarize the results in a dataframe
  results_df <- as.data.frame(do.call(rbind, all_results))
  # Give better names to the parametric quantile columns of the dataframe
  lower_quantile_par_names <- paste0("oneparQ", lower_q_levs, "geq0")
  upper_quantile_par_names <- paste0("oneparQ", upper_q_levs, "leq0")
  pinball_loss_par_names <- paste0("pinballOnepar",c(lower_q_levs, upper_q_levs))
  # Give better names to the nonparametric quantile columns of the dataframe
  lower_quantile_nonpar_names <- paste0("nonparQ", lower_q_levs, "geq0")
  upper_quantile_nonpar_names <- paste0("nonparQ", upper_q_levs, "leq0")
  pinball_loss_nonpar_names <- paste0("pinballNonpar", c(lower_q_levs, upper_q_levs))
  # Give better names to the linear quantile columns of the dataframe
  lower_quantile_lin_names <- paste0("linQ", lower_q_levs, "geq0")
  upper_quantile_lin_names <- paste0("linQ", upper_q_levs, "leq0")
  pinball_loss_lin_train_names <- paste0("pinballLin", c(lower_q_levs, upper_q_levs), "Train")
  pinball_loss_lin_test_names <- paste0("pinballLin", c(lower_q_levs, upper_q_levs), "Test")
  # mcqrnn columns
  lower_quantile_mcqrnn_names <- paste0("qrnn", lower_q_levs, "geq0")
  upper_quantile_mcqrnn_names <- paste0("qrnn", upper_q_levs,"leq0")
  pinball_loss_mcqrnn_train_names <- paste0("pinballQrnn", c(lower_q_levs, upper_q_levs), "Train")
  pinball_loss_mcqrnn_test_names <- paste0("pinballQrnn", c(lower_q_levs, upper_q_levs), "Train")
  colnames(results_df) <- c("num_samples",
                            "nu",
                            "dim",
                            "tau_min",
                            "tau_max",
                            "param_cond_func_idx",
                            lower_quantile_par_names,
                            upper_quantile_par_names,
                            pinball_loss_par_names,
                            lower_quantile_nonpar_names,
                            upper_quantile_nonpar_names,
                            pinball_loss_nonpar_names,
                            lower_quantile_lin_names,
                            upper_quantile_lin_names,
                            pinball_loss_lin_train_names,
                            pinball_loss_lin_test_names,
                            lower_quantile_mcqrnn_names,
                            upper_quantile_mcqrnn_names,
                            pinball_loss_mcqrnn_train_names,
                            pinball_loss_mcqrnn_test_names,
                            "train set loss",
                            "train set accuracy",
                            "test set loss",
                            "test set accuracy",
                            "MC_Integral",
                            "log_lik_simp_train_orig",
                            "log_lik_NN_train_orig",
                            "log_lik_true_train",
                            "log_lik_simp_test_orig",
                            "log_lik_NN_test_orig",
                            "log_lik_true_test",
                            "NN_num_params",
                            "simp_cop_num_params",
                            "deg_free",
                            "LRT_stat_train",
                            "LRT_stat_test",
                            "LRT_stat_overall",
                            "p_value_train",
                            "p_value_test",
                            "p_value_overall",
                            "AIC_NN",
                            "AIC_simp",
                            "BIC_NN",
                            "BIC_simp",
                            "wilcox_test_statistic",
                            "wilcox_p_value",
                            "wilcox_conf_int_bottom",
                            "wilcox_conf_int_top",
                            "t_test_statistic",
                            "t_test_p_value",
                            "t_test_conf_int_bottom",
                            "t_test_conf_int_top")
  if(filename != ""){
    write.csv(results_df, file = filename, row.names = FALSE)
    print(paste("Data saved to:", filename))
  }
  return(results_df)
}
# Run big simulation
# 2 different number of samples.
num_samples <- list(1000)
# fraction of number of noise to true samples. Higher, if only 1000 samples given.
nus <- list(10)
# dimensions of the simulated data.
dims <- list(3,5)
num_replications <- 10
struct_mats <- list(struct_mat_3d, struct_mat_5d)
families <- list(families_3d, families_5d)
initial_params <- list(params_3d, params_5d)
rotations <- list(rotations_3d, rotations_5d)

for(i in 1:num_replications){
  print(paste0("Replication ", i, " of ", num_replications))
  replication_name <- paste0("Repl",i)
  # Current data (Alternative: Use Sys.Date())
  current_date <- "20250623OnePar"
  # Construct the file name with the date
  num_samples_filename <- paste0(unlist(num_samples), collapse = ",")
  results_filename <- paste0("results/",
                             current_date,
                             "NSamples",
                             num_samples_filename,
                             replication_name,
                             ".csv")
  result_df <- run_simulations(
    num_samples=num_samples,
    dims=dims,
    nus=nus,
    tau_limits=tau_limits,
    struct_mats=struct_mats,
    families=families,
    initial_params=initial_params,
    rotations=rotations,
    filename=results_filename
  )
}


# # Plot Pairs copula data
# source("R/plotting_methods.R")
# #' plot pairs copula data
# #' @param num_samples: list of integers, that determine for what different
# #' number of samples drawn from a non-simplified vine copula the experiments should be run.
# #' @param dims: list of integers, that determine for what dimensions to run the tests
# #' (currently implemented for 3, 4 and 5)
# #' @param tau_limits: list of 2 dimensional vectors: For each of those vectors,
# #' the first entry determines, what the lower threshold of kendalls tau is,
# #' the second entry determines what the upper threshold of kendalls tau is in the
# #' conditional copulas
# #' @param struct_mats: List of regular vine matrices: Determine what structure to use.
# #' The entry in position i needs to have the same dimension as dims[[i]].
# #' @param families: List of (list of list of string): The i-th element
# #' contains the copula families corresponding to the copulas defined with struct_mats[[i]]
# #' @param initial_params: List of vectors: The i-th entry contains the parameters
# #' for the i-th unconditioned copula in the first tree, defined by struct_mats[[i]]
# #' @param rotations: List of (list of list of int): The i-th element contains the
# #' rotations of the copulas specified in struct_mats[[i]]
# #' @returns result_df: A Dataframe with the results of the experiments.
# plot_data <- function(num_samples,
#                             dims,
#                             tau_limits,
#                             struct_mats,
#                             families,
#                             initial_params,
#                             rotations
#                             ){
#   # Run big simulation
#   total_runs <- 0
#   # Initialize an empty list to store the results of each inner loop iteration
#   all_plots <- list()
#   for(sample_idx in 1:length(num_samples)){
#     # Check, if any of the dimensions are not implemented. If so stop execution and print.
#     incorrect_dimensions <- setdiff(unlist(dims), c(3,4,5))
#     if (length(incorrect_dimensions) > 0) {
#       print("The following dimensions are not implemented: ")
#       print(incorrect_dimensions)
#       print("Please remove these dimensions, then try to execute this function again.")
#       break
#     }
#     for(tau_idx in 1:length(tau_limits)){
#       tau_lower=tau_limits[[tau_idx]][1]
#       tau_upper=tau_limits[[tau_idx]][2]
#       for(dim_idx in 1:length(dims)){
#         # initialize param_cond_funcs (overwritten below)
#         param_cond_func_var <- -1
#         # Define the conditional parameter functions
#         if(dims[[dim_idx]] ==3){
#           param_cond_func_3d_lin <- list(
#             list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_func_3d_non_lin <- list(
#             list(u_to_param_quadratic(c(1,1),
#                                       tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_func_3d_non_lin_cubic <- list(
#             list(u_to_param_cubic(c(1,1.2,0.5),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_func_var <- list(
#             param_cond_func_3d_lin,
#             param_cond_func_3d_non_lin,
#             param_cond_func_3d_non_lin_cubic)
#         } else if(dims[[dim_idx]] ==4){
#           param_cond_funcs_4d_lin <- list(
#             list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_funcs_4d_non_lin <- list(
#             list(
#               u_to_param_quadratic(c(1,1.4),
#                                    tau_lower=tau_lower, tau_upper=tau_upper),
#               u_to_param_quadratic(c(1.2,0.9),
#                                    tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(
#               u_to_param_quadratic(c(-0.7,1.5,0.9,-1.2,1.5),
#                                    tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_funcs_4d_non_lin_cubic <- list(
#             list(u_to_param_cubic(c(0.7,1.2,-0.8),
#                                   tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_cubic(c(1.4,-1.8,1.0),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_cubic(c(0.7,1.8, -0.9,1.4,0.8,1.1,-1.2,-1.3,1.7),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_func_var <- list(
#             param_cond_funcs_4d_lin,
#             param_cond_funcs_4d_non_lin,
#             param_cond_funcs_4d_non_lin_cubic)
#         } else if(dims[[dim_idx]] ==5){
#           param_cond_funcs_5d_lin <- list(
#             list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_linear(c(0.2,0.5,0.3), tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_funcs_5d_non_lin <- list(
#             list(
#               u_to_param_quadratic(c(1,1),
#                                    tau_lower=tau_lower, tau_upper=tau_upper),
#               u_to_param_quadratic(c(1,0.7),
#                                    tau_lower=tau_lower, tau_upper=tau_upper),
#               u_to_param_quadratic(c(1,2),
#                                    tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_quadratic(c(0.7,0.5, 0.9,-1.2,1.5),
#                                       tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_quadratic(c(0.4,-0.8,1.0,0.7,2.0),
#                                       tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_quadratic(c(0.7,0.4,-0.9,1.3,0.75,1.3,-0.6,0.5,-1.1),
#                                       tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_funcs_5d_non_lin_cubic <- list(
#             list(u_to_param_cubic(c(1.0,1.0,0.5),
#                                   tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_cubic(c(1.0,0.7,0.3),
#                                   tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_cubic(c(1.0,2.0,1.1),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_cubic(c(0.7,0.8, -0.9,0.4,0.8,1.1,-1.2,-1.3,1.0),
#                                   tau_lower=tau_lower, tau_upper=tau_upper),
#                  u_to_param_cubic(c(0.4,0.6,1,0.7,2, -1.1,1.2,-0.9,-0.3),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)),
#             list(u_to_param_cubic(c(0.7, 0.5, -1.3,
#                                     1,   1.4, 1  ,-1.4,0.7,-0.8,
#                                     1.1,-1.4,0.7,-0.3,0.4,
#                                     0.8, -0.7,1.4,-1.2,-0.6),
#                                   tau_lower=tau_lower, tau_upper=tau_upper)))
#           param_cond_func_var <- list(
#             param_cond_funcs_5d_lin,
#             param_cond_funcs_5d_non_lin,
#             param_cond_funcs_5d_non_lin_cubic)
#         }
#         # Define the other necessary parameters for simulation
#         struct_mat_var <- struct_mats[[dim_idx]]
#         families_var <- families[[dim_idx]]
#         params_var <- initial_params[[dim_idx]]
#         rotations_var <- rotations[[dim_idx]]
#         for(par_idx in 1:length(param_cond_func_var)){
#           total_runs <- total_runs+1
#           print(paste0("Executing... ",
#                        total_runs,
#                        " of a total of ",
#                        length(dims)*length(num_samples)*length(tau_limits)*length(param_cond_func_var)))
#           non_simp_data <- simulate_non_simplified(n_samples = num_samples[[sample_idx]],
#                                                       struct = struct_mat_var,
#                                                       families=families_var,
#                                                       params = params_var,
#                                                       param_cond_funcs = param_cond_func_var[[par_idx]],
#                                                       rotations = rotations_var)
#           orig_data <- as.matrix(non_simp_data)
#           orig_data <- unname(orig_data)
#           param_cond_func_name <- "linear"
#           if(par_idx == 2){
#             param_cond_func_name <- "quadratic"
#           } else if (par_idx ==3){
#             param_cond_func_name <- "cubic"
#           }
#           all_plots[[length(all_plots)+1]] <- ggplotify::as.ggplot(
#             copula_pairs_ggplot(orig_data, title=paste0(param_cond_func_name, " tau_max= ", tau_upper))
#           )
#         }
#       }
#     }
#   }
#
#   return(all_plots)
# }
# # PAIRS COPULA PLOT
# n_samples_plot <- list(10000)
#
# # dimensions of the simulated data.
# dims <- list(3)
# # each list element contains first the lower, then the upper limit
# tau_limits <- list(c(0.001,0.3), c(0.001,0.6), c(0.001,0.9))
# # Define 3d parameters
# struct_mat_3d <- matrix(c(1,1,1,
#                           2,2,0,
#                           3,0,0), ncol=3, byrow=TRUE)
# families_3d <- list(list("frank", "gaussian"), list("frank"))
# params_3d <- list(c(ktau_to_par(family=families_3d[[1]][[1]], tau=0.2)),
#                   c(ktau_to_par(family=families_3d[[1]][[2]], tau=0.4)))
# rotations_3d <- list(list(0,0),0)
# # Define 5d parameters
# struct_mat_5d <- matrix(c(2,3,2,1,1,
#                           3,2,1,2,0,
#                           1,1,3,0,0,
#                           4,4,0,0,0,
#                           5,0,0,0,0), ncol=5, byrow=TRUE)
# families_5d <- list(list("frank", "gaussian","gaussian","frank"),
#                     list("frank","gaussian","gaussian"),
#                     list("gaussian", "frank"),
#                     list("gaussian"))
# params_5d <- list(c(ktau_to_par(family=families_5d[[1]][[1]], tau=0.2)),
#                   c(ktau_to_par(family=families_5d[[1]][[2]], tau=0.3)),
#                   c(ktau_to_par(family=families_5d[[1]][[3]], tau=0.4)),
#                   c(ktau_to_par(family=families_5d[[1]][[4]], tau=0.1)))
# rotations_5d <- list(list(0,0,0,0),list(0,0,0), list(0,0), list(0))
#
# struct_mats <- list(struct_mat_3d)
# families <- list(families_3d)
# initial_params <- list(params_3d)
# rotations <- list(rotations_3d)
# result_plots <- plot_data(num_samples = n_samples_plot,
#                           dims=dims,
#                           tau_limits=tau_limits,
#                           struct_mats = struct_mats,
#                           families=families,
#                           initial_params=initial_params,
#                           rotations=rotations)
# library(patchwork)
# length(result_plots)
# all_pairs_plots <- (result_plots[[1]] + result_plots[[2]] + result_plots[[3]]) /
#   (result_plots[[4]] + result_plots[[5]] + result_plots[[6]]) /
#   (result_plots[[7]] + result_plots[[8]] + result_plots[[9]])
# # all_pairs_plots
# ggplot2::ggsave(
#   filename = paste0("3d_non_simplified_all_configs_10000_plot.png"),
#   plot = all_pairs_plots,
#   width = 15,
#   height = 15,
#   dpi = 300,
#   bg="white"
# )
# # END PAIRS COPULA PLOT
# # Current data (Alternative: Use Sys.Date())
# current_date <- "20250616"
# # Construct the file name with the date
# num_samples_filename <- paste0(unlist(num_samples), collapse = ",")
# results_filename <- paste0("results/", current_date,"nSamples", num_samples_filename, ".csv")
# result_df <- run_simulations(
#   num_samples=num_samples,
#   dims=dims,
#   nus=nus,
#   tau_limits=tau_limits,
#   struct_mats=struct_mats,
#   families=families,
#   initial_params=initial_params,
#   rotations=rotations,
#   filename=results_filename
#   )
# head(result_df)
#
# q_levs_lower_for_latex <- c(0.05,0.10,0.15,0.2)
# q_levs_lower_latex_names <- paste0("q", q_levs_lower_for_latex, "geq0")
# q_levs_upper_for_latex <- c(0.95,0.9,0.85,0.8)
# q_levs_upper_latex_names <- paste0("q", q_levs_upper_for_latex, "geq0")
# cols_for_latex <-c("num_samples",
#                    "dim",
#                    "tau_lower",
#                    "tau_upper",
#                    "param_cond_func_idx",
#                    q_levs_lower_latex_names,
#                    q_levs_upper_latex_names,
#                    "train set loss",
#                    "test set loss",
#                    "AIC_NN",
#                    "AIC_simp",
#                    "log_lik_simp_train_orig",
#                    "log_lik_NN_train_orig",
#                    "log_lik_simp_test_orig",
#                    "log_lik_NN_test_orig",
#                    "NN_num_params",
#                    "simp_cop_num_params"
#                    )
# latex_result_df <- result_df[cols_for_latex]
# head(latex_result_df)
# latex_filename <- "results/LatexRelevantFields20250527.csv"
# write.csv(latex_result_df, file = latex_filename, row.names = FALSE)
# print(paste("Data saved to:", latex_filename))
# Small test
# result_df <- run_simulations(
#   num_samples=list(1000),
#   dims=list(3),
#   nus=list(10),
#   tau_limits=list(c(0.001,0.9)),
#   struct_mats=list(struct_mat_3d),
#   families=list(families_3d),
#   initial_params=list(params_3d),
#   rotations=rotations
#   )
#
# head(result_df)

# boxplot(result_df[,"q(0.05)>0"]/result_df[,"num_samples"])
# orig_data <- as.matrix(u_test)
# orig_data <- unname(orig_data)
# nu = 10
# fitted_vine <-vinecop(orig_data,family_set="onepar")
# # simulate from the simplified vine, to train a classifier later.
# num_samples <- nu * nrow(orig_data) # number of samples to create from the fitted vine
# simplified_samples <- rvinecop(num_samples, fitted_vine)
# split_output <- train_test_split(orig_data, simplified_samples, train_perc=0.8, stratified=TRUE)
# x_train <- split_output[[1]]
# x_test <- split_output[[2]]
# y_train <- split_output[[3]]
# y_test <- split_output[[4]]
# sum(y_train) /length(y_train)
# model <- 0
# model <- build_model(input_dim=ncol(x_train), use_tanh=TRUE)
# model %>% summary
# # initial evaluation (to test implementation)
# pre_training_eval <- model %>% evaluate(x_train, y_train)
# pre_training_eval
# model <- train_model(model, x_train, y_train, num_epochs=100)
# test_set_eval <- model %>% evaluate(x_test, y_test)
# test_set_eval[["loss"]]
# test_set_eval[["binary_accuracy"]]
# g_vals <- compute_gvals(model, orig_data, nu=nu)
# lower_q_levs <- seq(0.01,0.2,0.01)
# lower_q_levs
# upper_q_levs <- seq(0.8,0.99,0.01)
# upper_q_levs
# output <- perform_quant_reg(g_vals,
#                             orig_data,
#                             bottom_quantile_levels=lower_q_levs,
#                             top_quantile_levels = upper_q_levs)
# print(output)
# print(output[[1]] / nrow(orig_data))
# plot(lower_q_levs, output[[1]])
# plot(upper_q_levs, output[[2]])
