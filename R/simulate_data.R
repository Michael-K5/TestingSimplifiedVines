# This script simulates non-simplified copula data, using the functions from
# the script simulate_non_simplified_vine.R
source("R/simulate_non_simplified_vine.R")

# Simulate Data
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
tau_lower = -0.92
tau_upper = 0.92
param_cond_funcs_test <- list(list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_linear(c(0.2,0.5,0.3), tau_lower=tau_lower, tau_upper=tau_upper)))
u_test <- simulate_non_simp_parallel(n_samples = 50000,
                                  struct = struct_mat,
                                  families=family_test,
                                  params = params_test,
                                  param_cond_funcs = param_cond_funcs_test,
                                  rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
#head(u_test)
#fit.struct_mat<-vinecop(u_test,family_set="onepar",structure=struct_mat)
#print.data.frame(summary(fit.struct_mat),digit=2)
pairs_copula_data(u_test)
# plot_contours_1d(param_cond_func_1d=u_to_param_linear(c(1)), family_vals=c("frank", "gaussian", "joe"))
# plot_contours_2d(param_cond_func_2d=u_to_param_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper),
#                 family_name="gaussian")
# plot_contours_2d(param_cond_func_2d=u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper),
#                  family_name="gumbel")

# TEST as of 16.05.2025
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
pairs_copula_data(u_test)
# END TEST 16.05.2025

# TEST 2 as of 16.05.2025
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
tau_lower = -0.92
tau_upper = 0.92
param_cond_funcs_test <- list(list(u_to_param_non_lin_cub(c(1,1,0.5),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin_cub(c(1,0.7,0.3),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin_cub(c(1,2,0.1),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_non_lin_cub(c(0.7,0.3, 0.9,-0.4,0.8,0.1,0.2,0.3,0.4),
                                                      tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_non_lin_cub(c(0.4,0.6,1,0.7,2,0.01,0.02,0.04,0.03),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_non_lin_cub(c(0.2,0.5,0.3,
                                                            1,2,1,0.4,0.7,0.8,
                                                            0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1),
                                                      tau_lower=tau_lower, tau_upper=tau_upper)))
u_test <- simulate_non_simp_parallel(n_samples = 5000,
                                     struct = struct_mat,
                                     families=family_test,
                                     params = params_test,
                                     param_cond_funcs = param_cond_funcs_test,
                                     rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
pairs_copula_data(u_test)
# END TEST 2 of 16.05.2025

# Simulate different data
struct_mat_1 <- matrix(c(3,3,3,3,
                       2,2,2,0,
                       4,4,0,0,
                       1,0,0,0), ncol=4, byrow=TRUE)
families_1 <- list(list("frank", "gumbel","joe"),
                   list("clayton","gaussian"),
                   list("frank"))
params_1 <- list(c(ktau_to_par(family=families_1[[1]][[2]], tau=-0.2)),
                 c(ktau_to_par(family=families_1[[1]][[2]], tau=0.3)),
                 c(ktau_to_par(family=families_1[[1]][[3]], tau=0.1)))
param_cond_funcs_1 <- list(list(u_to_param_linear(c(1)), u_to_param_quadratic(c(1))),
                           list(u_to_param_linear(c(0.4,0.6))))
u_test_1 <- simulate_non_simp_parallel(n_samples = 8000,
                                     struct = struct_mat_1,
                                     families=families_1,
                                     params = params_1,
                                     param_cond_funcs = param_cond_funcs_1,
                                     rotations = list(list(0,0,0),list(0,0), list(0)))
#fit.struct_mat<-vinecop(u_test_1,family_set="onepar",structure=struct_mat)
#print.data.frame(summary(fit.struct_mat),digit=2)
pairs_copula_data(u_test_1)

struct_mat_2 <- matrix(c(1,1,1,
                         2,2,0,
                         3,0,0), ncol=3, byrow=TRUE)
u_test_2 <- simulate_non_simp_parallel(n_samples=4000,
                                    struct = struct_mat_2,
                                    families= list(list("frank", "frank"), list("frank")),
                                    params=list(c(1.3), c(2)),
                                    param_cond_funcs = list(list(u_to_param)), #u_to_param
                                    rotations=list(list(0,0), list(0)))
pairs_copula_data(u_test_2)
# see how a fit of a simplified copula looks like
#simplified_fit<-vinecop(u_test_2,family_set="onepar",structure=struct_mat_2)
# print.data.frame(summary(fit.struct_mat_2),digit=2)
#temp_test <- rvinecop(2000, simplified_fit)
#pairs_copula_data(temp_test)

# Simulate data of arbitrary dimension
data_dim <- 5
struct_mat_3 <- matrix(rep(0,data_dim^2), ncol=data_dim)
params_3 <- list()
families_3 <- list()
param_cond_funcs_3 <-list()
rotations_3 <- list()
param_vec <- runif(data_dim, min=1.1, max=3.0)
for(i in 1:data_dim){
  temp_fam <- list()
  temp_param_funcs <- list()
  temp_rotations <- list()
  params_3[[i]] <- c(param_vec[i])
  for(j in 1:(data_dim - i+1)){
    struct_mat_3[i,j] <- i
    temp_fam <- c(temp_fam, "frank")
    temp_param_funcs <- c(temp_param_funcs, u_to_param)
    temp_rotations <- c(temp_rotations, 0)
  }
  families_3[[i]] <- temp_fam
  param_cond_funcs_3[[i]] <- temp_param_funcs
  rotations_3[[i]] <- temp_rotations
}
u_test_3 <- simulate_non_simp_parallel(n_samples=5000,
                                    struct=struct_mat_3,
                                    families=families_3,
                                    params=params_3,
                                    param_cond_funcs=param_cond_funcs_3,
                                    rotations = rotations_3)
pairs_copula_data(u_test_3)
#fit.struct_mat<-vinecop(u_test_3,family_set="onepar",structure=struct_mat_3)
#print.data.frame(summary(fit.struct_mat),digit=2)

data_to_use <- u_test
# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()
data_dim <- paste0(ncol(data_to_use), "d_")
# Construct the file name with the date
csv_file_name <- paste0("data/non_simplified_sim_",data_dim, current_date, ".csv")
# Save as CSV
write.csv(data_to_use, file = csv_file_name, row.names = FALSE)

print(paste("Data saved to", csv_file_name, "\n"))

# Example how to read the data afterwards:
# temp <- read.csv(csv_file_name, header=TRUE)
# temp_mat <- as.matrix(temp)
# nrow(temp_mat)


