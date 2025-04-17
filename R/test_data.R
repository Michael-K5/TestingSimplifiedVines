source("R/simulate_non_simplified_vine.R")
# T values should be between T_min and T_max
T_min
T_max
fisher_z
# TEST
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
struct_mat <- matrix(c(3,3,3,3,
                       2,2,2,0,
                       4,4,0,0,
                       1,0,0,0), ncol=4, byrow=TRUE)
u_test <- simulate_non_simplified(n_samples = 3000,
                                  struct = struct_mat,
                                  families=list(list("frank", "frank","frank"), list("frank","frank"), list("frank")),
                                  params = list(c(2), c(1.3), c(1)),
                                  param_cond_func = u_to_param, #for tests: function(u, family) 3
                                  rotations = list(list(0,0,0),list(0,0), list(0)))
fit.struct_mat<-vinecop(u_test,family_set="onepar",structure=struct_mat)
print.data.frame(summary(fit.struct_mat),digit=2)
pairs_copula_data(u_test)

struct_mat_2 <- matrix(c(1,1,1,
                         2,2,0,
                         3,0,0), ncol=3, byrow=TRUE)
u_test_2 <- simulate_non_simplified(n_samples=5000,
                                    struct = struct_mat_2,
                                    families= list(list("frank", "frank"), list("frank")),
                                    params=list(c(1.3), c(2)),
                                    param_cond_func = u_to_param,
                                    rotations=list(list(0,0), list(0)))
pairs_copula_data(u_test_2)
# see how a fit of a simplified copula looks like
#simplified_fit<-vinecop(u_test_2,family_set="onepar",structure=struct_mat_2)
# print.data.frame(summary(fit.struct_mat_2),digit=2)
#temp_test <- rvinecop(2000, simplified_fit)
#pairs_copula_data(temp_test)

# Get the current date in YYYY-MM-DD format
current_date <- Sys.Date()

# Construct the file name with the date
csv_file_name <- paste0("data/non_simplified_sim_", current_date, ".csv")

# Save as CSV
write.csv(u_test_2, file = csv_file_name, row.names = FALSE)

print(paste("Data saved to", csv_file_name, "\n"))

# Example how to read the data afterwards:
# temp <- read.csv(csv_file_name, header=TRUE)
# temp_mat <- as.matrix(temp)
# nrow(temp_mat)
