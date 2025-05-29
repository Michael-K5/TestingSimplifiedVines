# Script to run tests
# TEST
source("R/simulate_non_simplified_vine.R")
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

plot_contours_1d()
plot_contours_2d()
plot_contours_2d(family_name="frank", param_cond_func_2d=u_to_param_linear(c(0.1,0.9)), manual_subtitle="Weights of linear function: 0.1,0.9")
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
pairs_copula_data_custom(u_test)
pairs_copula_data_custom(u_test, max_samples=2000, plot_emp_ktau=TRUE)

# own implementation of pairs_copula_data
df <- data.frame(u_test)

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

temp_struct <- matrix(c(3, 2, 3 , 4 , 4 ,
2 , 3 , 4 , 3 , 0 ,
4 , 4 , 2 , 0 , 0 ,
1 , 1 , 0 , 0 , 0 ,
5 , 0 , 0 , 0 , 0), ncol=5, byrow=TRUE)
plot(rvine_matrix(temp_struct),1:4)


# Test tanh
x <- -10000:10000 / 10
y1 <- tanh(x)
y2 <- inverse_fisher_transform(x)
plot(y1,y2)
par(mfrow=c(1,1))
plot(x,y1)
plot(x,y2)
plot(x,atanh(y1))

sigmoid_fun <- function(temp){
  return(1/(1+exp(-temp)))
}

sigmoid_inv <- function(temp){
  return(-log(1/temp  -1))
}

x <- -1000:1000 / 10
plot(x, sigmoid_inv(sigmoid_fun(x)))
sigmoid_fun(10)
sigmoid_inv(sigmoid_fun(10))
sigmoid_fun(10000)

y <- c(0.1,0.2,0.6)
qnorm(y)

matrix(poly(y, degree=2,simple=TRUE,raw=TRUE))
a <- 2
b <- 3
data <- data.frame(a = a, b = b)
model.matrix(~ poly(a, b, degree = 2, raw = TRUE), data)

# Test plot the kendalls tau values for different conditioning values

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1).
#' That function takes the dot product between a and u, scales that to
#' T_values and calculates a kendall's tau value using the tanh.
#' @param a: a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
#' @param tau_lower: Lowest tau value, defaults to -0.92
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value, defaults to 0.92
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u, which calculates the kendalls tau value of
#' a copula given the conditioned values.
u_to_ktau_linear <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u){
    tryCatch({
      T_upper <- fisher_z_transform(tau_upper)
      T_lower <- fisher_z_transform(tau_lower)
      T_val <- T_lower
      if(length(u) == 1){
        arg <- u
        T_val <- (T_upper- T_lower) * arg + T_lower
      } else {
        arg <- a %*% u
        T_val <- (T_upper - T_lower) * arg + T_lower
      }
      tau <- inverse_fisher_transform(T_val)
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
    })
  })
}

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function),
#' builds polynomial terms of degree 2 of these terms and then transforms them
#' to the interval (tau_lower, tau_upper) using a scaled tanh.
#' If dim(u) <=3 includes all terms u[i]*u[j], (powers and mixed terms).
#' If dim(u) > 3 only includes the individual powers, u[i]^k, k=1,2.
#' @param a: a vector of weights
#' @param tau_lower: Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_ktau_quadratic <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      # define parameters for the scaled tanh, so that the result will be between
      # tau_lower and tau_upper
      scaling_factor <- (tau_upper - tau_lower) / 2
      shift <- (tau_upper + tau_lower) / 2
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        scaling_factor <- (tau_upper - 0.001) / 2
        shift <- (tau_upper + 0.001) / 2
      }
      arg <- qnorm(u)
      tau <- 0
      if(length(u) == 1){
        temp <- c(arg, arg^2) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==2) {
        temp <- c(arg[1], arg[2], arg[1]^2, arg[2]^2, arg[1]*arg[2]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==3){
        temp <- c(arg[1], arg[2], arg[3],
                  arg[1]^2, arg[2]^2, arg[3]^2,
                  arg[1]*arg[2], arg[1]*arg[3], arg[2]*arg[3]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      }else{
        # do not include mixed terms here.
        arg <- poly(arg, 2, raw=TRUE) # raw to just evaluate the polynomial terms
        # only keep the entries and list them in a vector with
        # c(arg[1],arg[2],..., arg[1]^2, arg[2]^2,...)
        arg <- c(unname(arg[1:nrow(arg),]))
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
    })
  })
}

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function),
#' builds polynomial terms of degree 2 of these terms and then transforms them
#' to the interval (tau_lower, tau_upper) using a scaled tanh.
#' If dim(u) <=3 includes powers and mixed terms,
#' If dim(u) > 3 only includes the individual powers, u[i]^k, k=1,2,3.
#' @param a: a vector of weights
#' @param tau_lower: Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_ktau_cubic <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      # define parameters for the scaled tanh, so that the result will be between
      # tau_lower and tau_upper
      scaling_factor <- (tau_upper - tau_lower) / 2
      shift <- (tau_upper + tau_lower) / 2
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        scaling_factor <- (tau_upper - 0.001) / 2
        shift <- (tau_upper + 0.001) / 2
      }
      arg <- qnorm(u)
      tau <- 0
      if(length(u) == 1){
        temp <- c(arg, arg^2, arg^3) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==2) {
        temp <- c(arg[1], arg[2],
                  arg[1]^2, arg[2]^2, arg[1]*arg[2],
                  arg[1]^3, arg[2]^3,arg[1]^2 * arg[2], arg[1]*arg[2]^2) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==3){
        temp <- c(arg[1], arg[2], arg[3],
                  arg[1]^2, arg[2]^2, arg[3]^2,
                  arg[1]*arg[2], arg[1]*arg[3], arg[2]*arg[3],
                  arg[1]^3, arg[2]^3, arg[3]^3,
                  arg[1]^2*arg[2], arg[1]^2*arg[3],
                  arg[2]^2*arg[1], arg[2]^2*arg[3],
                  arg[3]^2*arg[1], arg[3]^2*arg[2], arg[1]*arg[2]*arg[3]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      }else{
        # do not include mixed terms here.
        arg <- poly(arg, 3, raw=TRUE) # raw to just evaluate the polynomial terms
        # only keep the entries and list them in a vector with
        # c(arg[1],arg[2],..., arg[1]^2, arg[2]^2,..., arg[1]^3,arg[2]^3,...)
        arg <- c(unname(arg[1:nrow(arg),]))
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
    })
  })
}

tau_lower = 0.001
tau_upper = 0.9
func_lin_1d <- u_to_ktau_linear(c(1), tau_lower=tau_lower, tau_upper=tau_upper)
func_lin_2d <- u_to_ktau_linear(c(0.7,0.3), tau_lower=tau_lower, tau_upper=tau_upper)
func_lin_3d <- u_to_ktau_linear(c(0.2,0.5,0.3), tau_lower=tau_lower, tau_upper=tau_upper)
func_quadratic_1d <- u_to_ktau_quadratic(c(1,1),tau_lower=tau_lower, tau_upper=tau_upper)
func_quadratic_2d <- u_to_ktau_quadratic(c(0.7,0.5, 0.9,-1.2,1.5),tau_lower=tau_lower, tau_upper=tau_upper)
func_quadratic_3d <- u_to_ktau_quadratic(
  c(0.7,0.4,-0.9,1.3,0.75,1.3,-0.6,0.5,-1.1),
  tau_lower=tau_lower, tau_upper=tau_upper)
func_cubic_1d <- u_to_ktau_cubic(c(1.0,1.0,0.5),
                                tau_lower=tau_lower, tau_upper=tau_upper)
func_cubic_2d <- u_to_ktau_cubic(c(0.7,0.8, -0.9,0.4,0.8,1.1,-1.2,-1.3,1.0),
                                tau_lower=tau_lower, tau_upper=tau_upper)
func_cubic_3d <- u_to_ktau_cubic(
  c(0.7, 0.5, -1.3,
    1,   1.4, 1  ,-1.4,0.7,-0.8,
    1.1,-1.4,0.7,-0.3,0.4,
    0.8, -0.7,1.4,-1.2,-0.6),
  tau_lower=tau_lower, tau_upper=tau_upper
)

#' plot of the cond to ktau function with 1 conditioning variable
plot_cond_to_ktau_1d <- function(func_list, titles=-1, u_cond_vals=1:99/100){
  par(mfrow=c(length(func_list), 1),
      mar = c(3, 3, 1.5, 1.5),
      oma = c(1.5, 0, 1.5, 0),
      mgp=c(2,1,0))
  for(i in 1:length(func_list)){
    ktau_vals <- sapply(u_cond_vals, func_list[[i]])
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    plot(u_cond_vals, ktau_vals, type="l", main=plot_title, xlab="u",ylab="ktau",cex.lab=0.9)
  }
  par(mfrow=c(1,1))
}
plot_cond_to_ktau_1d(list(func_lin_1d, func_quadratic_1d, func_cubic_1d),
                     titles = c("Linear function 1d",
                                "Quadratic function 1d",
                                "Cubic function 1d"))

#' 3d plot of the cond to ktau function, for 2 conditioning variables
plot_cond_to_ktau_2d <- function(func_list,
                                 titles=-1,
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100){
  for(i in 1:length(func_list)){
    ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)), ncol=length(u_cond_vals_2))
    for(j in 1:length(u_cond_vals_1)){
      for(k in 1:length(u_cond_vals_2)){
        ktau_vals[j,k] <- func_list[[i]](c(u_cond_vals_1[j], u_cond_vals_2[k]))
      }
    }
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    # Create a 3D surface plot
    plot3D::persp3D(u_cond_vals_1, u_cond_vals_2, ktau_vals,
                    theta = 30, phi = 20, axes=TRUE, ticktype="detailed",
                    xlab="u1", ylab="u2", zlab="ktau",
                    col = "lightblue", border = "black", main=plot_title)
  }
}
plot_cond_to_ktau_2d(list(func_lin_2d),
                     titles=c("Linear function 2d"))
plot_cond_to_ktau_2d(list(func_quadratic_2d),
                     titles=c("Quadratic function 2d"))
plot_cond_to_ktau_2d(list(func_cubic_2d),
                     titles=c("Cubic function 2d"))

#' contours of the cond to ktau function with 2 conditional variables
plot_ktau_contours_2d <- function(func_list,
                             titles=-1,
                             u_cond_vals_1=1:99/100,
                             u_cond_vals_2=1:99/100){
  for(i in 1:length(func_list)){
    ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                        ncol=length(u_cond_vals_2))
    for(j in 1:length(u_cond_vals_1)){
      for(k in 1:length(u_cond_vals_2)){
        ktau_vals[j,k] <- func_list[[i]](c(u_cond_vals_1[j], u_cond_vals_2[k]))
      }
    }
    #ktau_vals <- outer(u_cond_vals_1, u_cond_vals_2, function(u1,u2) func_2d(c(u1,u2)))
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    # Create a 3D surface plot
    contour(u_cond_vals_1, u_cond_vals_2, ktau_vals, main=plot_title,
            xlab="u1", ylab="u2")
  }
}
plot_ktau_contours_2d(c(func_lin_2d, func_non_lin_2d, func_non_lin_cubic_2d),
                      titles=c("Linear function 2d",
                               "Non-linear function 2d",
                               "Non-linear cubic function 2d"))


#' Several 3d Plots for in total 3 conditioning variables, always one is fixed for the plots.
plot_cond_to_ktau_3d <- function(cond_to_ktau_func,
                                 title="",
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100,
                                 u_cond_vals_fixed=c(0.25,0.5,0.75)){
  # order of mar and oma: bottom, left, top, right. oma=outer margin, mar = inner margin
  par(mfrow = c(3, length(u_cond_vals_fixed)), mar = c(0.5, 0.5, 1, 0.5), oma = c(0.5, 0.5, 4, 0.5))
  for(fixed_dim in 1:3){
    for(i in 1:length(u_cond_vals_fixed)){
      ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                          ncol=length(u_cond_vals_2))
      for(j in 1:length(u_cond_vals_1)){
        for(k in 1:length(u_cond_vals_2)){
          if(fixed_dim==1){
            # dimension 1 is fixed
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_fixed[i], u_cond_vals_1[j], u_cond_vals_2[k]))
          } else if(fixed_dim==2){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_fixed[i], u_cond_vals_2[k]))
          } else if(fixed_dim==3){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_2[k], u_cond_vals_fixed[i]))
          }
        }
      }
      # Determine title and axis labels
      plot_title <- paste0(
        "u", fixed_dim,"=", u_cond_vals_fixed[i])
      xlab_text <- "u1"
      ylab_text <- "u2"
      if(fixed_dim == 1){
        xlab_text <- "u2"
        ylab_text <- "u3"
      } else if(fixed_dim==2){
        xlab_text <- "u1"
        ylab_text <- "u3"
      }
      # Create a 3d plot. (theta and phi determine the viewing angles)
      plot3D::persp3D(u_cond_vals_1, u_cond_vals_2, ktau_vals,
                      theta = 30, phi = 20, axes=TRUE, ticktype="detailed",
                      cex.axis=0.5, cex.lab =0.7,
                      xlab=xlab_text, ylab=ylab_text,zlab="ktau",
                      col = "lightblue", border = "black", main=plot_title)
    }
  }
  mtext(title, outer=TRUE, cex=1, line=2, font=2)
}
plot_cond_to_ktau_3d(func_lin_3d,
                     title="Linear Function 3d")
plot_cond_to_ktau_3d(func_quadratic_3d,
                     title="Quadratic function 3d")
plot_cond_to_ktau_3d(func_cubic_3d,
                     title="Cubic function 3d")

#' Contour plot of ktau for 3 conditioning variables, always one fixed
#' Several 3d Plots for in total 3 conditioning variables, always one is fixed for the plots.
plot_ktau_contour_3d <- function(cond_to_ktau_func,
                                 title="",
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100,
                                 u_cond_vals_fixed=c(0.25,0.5,0.7)){
  par(mfrow = c(3, length(u_cond_vals_fixed)), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(1.5, 4, 4, 4))
  for(fixed_dim in 1:3){
    for(i in 1:length(u_cond_vals_fixed)){
      ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                          ncol=length(u_cond_vals_2))
      for(j in 1:length(u_cond_vals_1)){
        for(k in 1:length(u_cond_vals_2)){
          if(fixed_dim==1){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_fixed[i], u_cond_vals_1[j], u_cond_vals_2[k]))
          } else if(fixed_dim==2){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_fixed[i], u_cond_vals_2[k]))
          } else if(fixed_dim==3){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_2[k], u_cond_vals_fixed[i]))
          }
        }
      }
      # Adjust title and axis labels
      plot_title <- paste0(
        "u", fixed_dim,"=", u_cond_vals_fixed[i])
      xlab_text <- "u1"
      ylab_text <- "u2"
      if(fixed_dim == 1){
        xlab_text <- "u2"
        ylab_text <- "u3"
      } else if(fixed_dim==2){
        xlab_text <- "u1"
        ylab_text <- "u3"
      }
      # Draw a contour plot
      contour(u_cond_vals_1, u_cond_vals_2, ktau_vals,
              xlab=xlab_text, ylab=ylab_text, main=plot_title)
    }
  }
  mtext(title, outer=TRUE, cex=1, line=2, font=2)
}
plot_ktau_contour_3d(func_lin_3d,
                     title="Linear Function 3d")
plot_ktau_contour_3d(func_non_lin_3d,
                     title="Non-linear function 3d",
                     u_cond_vals_fixed=c(0.25,0.5,0.75))
plot_ktau_contour_3d(func_non_lin_cubic_3d,
                     title="Non-linear cubic function 3d",
                     u_cond_vals_fixed=c(0.25,0.5,0.75))



# test for poly
args_list <- list(1, c(1,2), c(1,2,3))

cub <- poly(args_list[[2]], 3, raw=TRUE)
arg_to_mult <- c(unname(cub[1:nrow(cub),]))
arg_to_mult
print(arg_to_mult)
print(unname(cub[1:nrow(cub),]))



library(dplyr)
library(knitr)
latex_table <- read.csv("results/LatexRelevantFields20250527.csv")
# Assuming your data frame is called `df`
lower_q_levs_latex_table <- c(0.05,0.1,0.15,0.2)
upper_q_levs_latex_table <- c(0.95,0.9,0.85,0.8)
q_levs_latex_table <- c(lower_q_levs_latex_table, upper_q_levs_latex_table)
col_names_q_latex_table <- paste0("q.", q_levs_latex_table, "..0")

latex_result_df <- latex_table %>%
  mutate(across(all_of(col_names_q_latex_table), ~ .x / num_samples))

latex_result_df <- latex_result_df %>%
  mutate(across(all_of(col_names_q_latex_table), ~ round(.x, 4)))
latex_result_df <-  latex_result_df %>%
  mutate(param_cond_func_idx = recode(param_cond_func_idx, `1` = "linear", `2` = "quadratic", `3` = "cubic"))

latex_result_df <- latex_result_df %>%
  select(-all_of(c("train.set.loss", "test.set.loss")))
df_1000 <- latex_result_df %>%
  filter(abs(num_samples - 1000) < 1)
df_1000  <- df_1000  %>%
  select(-all_of(c("num_samples", "tau_lower")))
df_10000 <- latex_result_df %>%
  filter(abs(num_samples - 10000) < 1)
df_10000  <- df_10000  %>%
  select(-all_of(c("num_samples", "tau_lower")))

kable(df_1000, format="latex", booktabs="TRUE")
kable(df_10000,format="latex", booktabs="TRUE")
head(latex_result_df)
kable(latex_result_df, format = "latex", booktabs = TRUE)
# library(xtable)
# sink("results/latexTable20250528.tex")
# print(xtable(latex_result_df), include.rownames = FALSE, booktabs = TRUE)
# sink()

library(rvinecopulib)
struct_mat_4d <- matrix(c(3,3,3,3,
                          2,2,2,0,
                          4,4,0,0,
                          1,0,0,0), ncol=4, byrow=TRUE)

struct_mat_3d <- matrix(c(1,1,1,
                          2,2,0,
                          3,0,0), ncol=3, byrow=TRUE)
plot(as_rvine_matrix(struct_mat_3d),1:2)
plot(as_rvine_matrix(struct_mat_4d),1:3)
