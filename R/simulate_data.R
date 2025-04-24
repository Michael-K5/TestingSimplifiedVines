# This script simulates non-simplified copula data, using the functions from
# the script simulate_non_simplified_vine.R
source("R/simulate_non_simplified_vine.R")
# T values should be between T_min and T_max
T_min
#u_to_param(0.1, family="frank")
#u_to_param(0.7, family="frank")
T_max
# TEST
# test_tau <- 1:99 / 50 -1
# test_T <- fisher_z_transform(test_tau)
# test_tau_1 = inverse_fisher_transform(test_T)
# plot(test_tau, test_T)
# test_mat <- matrix(c(1,1,1,
#                      3,3,0,
#                      2,0,0)
#                    ,byrow=TRUE, ncol=3)
# test_result <- permute_indices(test_mat)
# test_result
# Test get_max_matrix
# temp <- matrix(c(1,1,1,1,
#                  3,2,2,0,
#                  2,3,0,0,
#                  4,0,0,0), nrow=4, byrow=TRUE)
# max_mat_test <- get_max_matrix(temp)
# print(max_mat_test)

# Simulate Data
struct_mat <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
u_test <- simulate_non_simp_parallel(n_samples = 5000,
                                  struct = struct_mat,
                                  families=list(list("frank", "frank","frank","frank"), list("frank","frank","frank"), list("frank", "frank"), list("frank")),
                                  params = list(c(2), c(1.3), c(1), c(1.5)),
                                  param_cond_funcs = list(list(u_to_param, u_to_param, u_to_param), list(u_to_param, u_to_param), list(u_to_param)), #for tests: function(u, family) 3
                                  rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
fit.struct_mat<-vinecop(u_test,family_set="onepar",structure=struct_mat)
print.data.frame(summary(fit.struct_mat),digit=2)
pairs_copula_data(u_test)

struct_mat_1 <- matrix(c(3,3,3,3,
                       2,2,2,0,
                       4,4,0,0,
                       1,0,0,0), ncol=4, byrow=TRUE)
u_test_1 <- simulate_non_simp_parallel(n_samples = 5000,
                                     struct = struct_mat_1,
                                     families=list(list("frank", "frank","frank"), list("frank","frank"), list("frank")),
                                     params = list(c(2), c(1.3), c(1)),
                                     param_cond_funcs = list(list(u_to_param, u_to_param), list(u_to_param)), #for tests: function(u, family) 3
                                     rotations = list(list(0,0,0),list(0,0), list(0)))
fit.struct_mat<-vinecop(u_test_1,family_set="onepar",structure=struct_mat)
print.data.frame(summary(fit.struct_mat),digit=2)
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


