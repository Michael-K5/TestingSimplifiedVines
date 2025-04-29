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
