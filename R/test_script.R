# Script to run tests
# TEST
# source("R/simulate_non_simplified_vine.R")
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

library(vinereg)
# define data
# testdata -> delete later.
u_obs <- data.frame(matrix(runif(400), ncol=4))
r_vals <- rep(0,100)
for (i in 1:nrow(u_obs)){
  r_vals[i] <- sum(u_obs[i,])
}
fit <- vinereg(r_vals ~ ., family="onepar",data=u_obs)
summary(fit)

plot_effects(fit)

quantiles <- predict(fit, u_obs, alpha=c(0.05,0.95))
for (i in 1:length(r_vals)){
  if(r_vals[i] < quantiles[i,1] | r_vals[i] > quantiles[i,2]){
    print(paste(i, "r_val outside 90% confidence interval"))
    print(r_vals[i])
    print(quantiles[i,])
  }
}


# test uniform margins
last_data_simulation_date <- "2025-04-23"
data_dim <- "4"
# load data
csv_filename <- paste0("data/non_simplified_sim_",data_dim,"d_",last_data_simulation_date,".csv")
orig_data <- as.matrix(read.csv(csv_filename))
orig_data <- unname(orig_data) #remove col- and rownames

for(i in 1:data_dim){
  hist(orig_data[,i], main = paste("Plot Number", i), probability=TRUE, breaks=seq(0,1,by=0.05))
}
breaks <- seq(0,1,by=0.1)
for(i in 1:data_dim){
  print(paste(i, "of", data_dim))
  observed <- table(cut(orig_data[,i], breaks = breaks, include.lowest = TRUE))

  # expected frequency under uniform distribution is same count in each bin
  expected <- rep(length(orig_data[,i]) / length(observed), length(observed))

  # chi-square test
  print(chisq.test(x = observed, p = rep(1/length(observed), length(observed))))
}

# quantile regression directly on the p-values -> not very good, as the quantiles are
# then sometimes not in the interval from 0 to 1.
obs_data <- data.frame(orig_data)
q_reg_p <- vinereg(predictions ~.,family_set="parametric", data=obs_data)
p_quantiles <- predict(q_reg_p, obs_data, alpha=c(0.1,0.9))
simp_not_sufficient <- 0
simp_sufficient <- 0
for( i in 1:nrow(p_quantiles)){
  if(p_quantiles[i,1] > 0.5){
    simp_not_sufficient <- simp_not_sufficient + 1
  } else if (p_quantiles[i,2] < 0.5){
    simp_sufficient <- simp_sufficient +1
  }
}
simp_not_sufficient # model can confidently assess that a sample is non-simplified
simp_sufficient # model can not detect that a sample is non-simplified
nrow(p_quantiles) -(simp_not_sufficient + simp_sufficient) # non conclusive observations
simp_not_sufficient / nrow(p_quantiles) #share of samples where simplified is not sufficient



# TEMP
source("R/simulate_non_simplified_vine.R")
library(rvinecopulib)
library(MASS)

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
params_test <- list(c(ktau_to_par(family=family_test[[1]][[1]], tau=-0.2)),
                    c(ktau_to_par(family=family_test[[1]][[2]], tau=0.3)),
                    c(ktau_to_par(family=family_test[[1]][[3]], tau=-0.1)),
                    c(ktau_to_par(family=family_test[[1]][[4]], tau=0.1)))
tau_lower = -0.9
tau_upper = 0.9
param_cond_funcs_test <- list(list(u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper),
                                   u_to_param_linear(c(0.4,0.6), tau_lower=tau_lower, tau_upper=tau_upper)),
                              list(u_to_param_linear(c(0.2,0.5,0.3), tau_lower=tau_lower, tau_upper=tau_upper)))
u_test <- simulate_non_simp_parallel(n_samples = 10000,
                                     struct = struct_mat,
                                     families=family_test,
                                     params = params_test,
                                     param_cond_funcs = param_cond_funcs_test,
                                     rotations = list(list(0,0,0,0),list(0,0,0), list(0,0), list(0)))
pairs_copula_data(u_test)

# own implementation of pairs_copula_data
df <- data.frame(u_test)

# Start the plotting
pairs_custom <- function(df) {
  n <- ncol(df)
  par(mfrow = c(n, n), mar = c(0, 0, 0, 0), oma = c(1.5, 1.5, 1.5, 1.5))

  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        hist(df[[i]], main = "", xlab = "", ylab = "",
             col = "darkgrey", border = "white", axes=FALSE)
        box()
        mtext(paste0("u",i), side = 3, line = -1, adj = 0.5, cex = 0.8, font=2)
      } else if (i < j) {
        plot(df[[j]], df[[i]], xlab = "", ylab = "", pch = ".", col = "black", axes=FALSE)
        box()
      } else {
        # Contour plot using 2D density estimation
        x <- df[[j]]
        y <- df[[i]]
        # convert to R using normal quantile function (yields easier to read contours)
        x <- qnorm(x)
        y <- qnorm(y)
        kde <- MASS::kde2d(x, y, n = 30, lims = c(-3, 3, -3, 3)) #lims = c(0, 1, 0, 1)
        contour(kde, nlevels=8,drawlabels = FALSE, xlab = "", ylab = "", axes = FALSE)
        box()
      }
    }
  }
}

# Run the custom plotting function
pairs_custom(df)


# Plot trees of copula structures
library(ggraph)
# M = matrix(c(1,7,6,7,7,7,7,
#              7,6,7,2,2,2,0,
#              2,2,2,6,6,0,0,
#              6,1,1,1,0,0,0,
#              5,5,5,0,0,0,0,
#              4,4,0,0,0,0,0,
#              3,0,0,0,0,0,0),ncol=7,nrow=7,byrow=TRUE)
# print(M)
#
# rvine_tree_struct= rvine_matrix(M)
# plot(rvine_tree_struct,1:6)# second argument tells R which trees to plot.


struct_mat <- matrix(c(2,3,2,1,1,
                       3,2,1,2,0,
                       1,1,3,0,0,
                       4,4,0,0,0,
                       5,0,0,0,0), ncol=5, byrow=TRUE)
rvine_tree_struct <- rvine_matrix(struct_mat)
plot(rvine_tree_struct,1:4)

d_vine_struct_5d <- matrix(c(4,3,2,1,1,
                             3,2,1,2,0,
                             2,1,3,0,0,
                             1,4,0,0,0,
                             5,0,0,0,0), ncol=5, byrow=TRUE)
dvine_tree_struct <- rvine_matrix(d_vine_struct_5d)
plot(dvine_tree_struct,1:4)
