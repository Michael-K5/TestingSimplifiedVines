# This script simulates non-simplified copula data, using the functions from
# the script simulate_non_simplified_vine.R
source("R/simulate_non_simplified_vine.R")
# T values should be between T_min and T_max
T_min
#u_to_param(0.1, family="frank")
#u_to_param(0.7, family="frank")
T_max

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function takes the dot product between a and u
#' and returns a linear function of that dot product.
#' @param a: a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
u_to_param_linear <- function(a){
  return(function(u, family="gaussian"){
    tryCatch({
      T_upper <- fisher_z_transform(tau_max)
      T_lower <- fisher_z_transform(tau_min)
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe")){
        T_lower <- fisher_z_transform(0.001)
      }
      T_val <- T_min
      if(length(u) == 1){
        arg <- u
        T_val <- (T_upper- T_lower) * arg + T_lower
      } else {
        arg <- a %*% u
        T_val <- (T_upper - T_lower) * arg + T_lower
      }
      tau <- inverse_fisher_transform(T_val)
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          The function u_to_param_linear only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function takes the dot product between a and u
#' and returns a quadratic function of that dot product.
#' @param a: a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
u_to_param_quadratic <- function(a){
  return(function(u, family="gaussian"){
    tryCatch({
      T_upper <- fisher_z_transform(tau_max)
      T_lower <- fisher_z_transform(tau_min)
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe")){
        T_lower <- fisher_z_transform(0.01)
      }
      T_val <- T_min
      if(length(u) == 1){
        arg <- u
        T_val <- 4*(T_upper- T_lower) * ((arg-0.5)^2) + T_lower
      } else {
        arg <- a %*% u
        T_val <- (T_upper - T_lower) * ((arg-0.5)^2) + T_lower
      }
      tau <- inverse_fisher_transform(T_val)
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          The function u_to_param_quadratic only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}


# Simulate Data
struct_mat <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
family_test <- list(list("frank", "clayton","gaussian","frank"),
                    list("frank","gaussian","joe"),
                    list("gaussian", "gumbel"),
                    list("gaussian"))
params_test <- list(c(ktau_to_par(family=family_test[[1]][[2]], tau=0.2)),
                    c(ktau_to_par(family=family_test[[1]][[2]], tau=-0.3)),
                    c(ktau_to_par(family=family_test[[1]][[3]], tau=0.1)),
                    c(ktau_to_par(family=family_test[[1]][[4]], tau=-0.1)))
param_cond_funcs_test <- list(list(u_to_param_linear(c(1)), u_to_param_linear(c(1)), u_to_param_linear(c(1))),
                              list(u_to_param_linear(c(0.7,0.3)), u_to_param_linear(c(0.4,0.6))),
                              list(u_to_param_linear(c(0.2,0.5,0.3))))
u_test <- simulate_non_simp_parallel(n_samples = 4000,
                                  struct = struct_mat,
                                  families=family_test,
                                  params = params_test,
                                  param_cond_funcs = param_cond_funcs_test,
                                  rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
#head(u_test)
#fit.struct_mat<-vinecop(u_test,family_set="onepar",structure=struct_mat)
#print.data.frame(summary(fit.struct_mat),digit=2)
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


