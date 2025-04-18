# This script contains functions to simulate from a Non-Simplified Vine Copula
library(rvinecopulib)
#' Computes the fisher Z transform (artanh)
#' @param tau: the value of Kendall's tau to transform to a real number.
#' Can also be a vector of taus.
#'
#' @return The fisher transformation(s) of tau. If tau is a number, this returns a number
#' if tau was a vector this returns a vector
fisher_z_transform = function(tau){
  return(1/2 * log((1+tau)/(1-tau)))
}
#' Computes the inverse of the fisher Z transform, which is the tanh
#' @param T_val: The T-value that is supposed to be turned into a Kendall's tau again
#'
#' @return The inverse fisher transformation(s) of T_val. If T is a number, this returns a number
#' if T is a vector this returns a vector
inverse_fisher_transform <- function(T_val){
  # computes the inverse of the fisher Z transform (tanh)
  return((exp(2*T_val) - 1)/(exp(2*T_val) + 1))
}

# tau should be in the range -0.97, +0.97 to avoid numerical problems later on
tau_min <- - 0.92
tau_max <- 0.92
T_min <- fisher_z_transform(tau_min)
T_max <- fisher_z_transform(tau_max)
#T_min
#T_max

#' maps a u value to the specific parameter for the third copula, by
#' finding a value on the fisher-z-transform-scale and mapping that back to a parameter
#' @param u: the samples u value (between 0 and 1)
#' @param family: the family of the third copula in the sample
u_to_param <- function(u, family="gaussian"){
  tryCatch({
    #arg = ifelse(length(u)==1,u,max(u))
    arg=ifelse(length(u)==1,u,mean(u))
    #T_val <- (T_max - T_min)*arg + T_min
    T_val <- ifelse(arg > 0.5,T_max,T_min)
    # T_val <- 4*(T_max - T_min) * ((arg-0.5)^2) + T_min
    tau <- inverse_fisher_transform(T_val)
    param <- ktau_to_par(family=family, tau=tau)
    return(param)
  },
  error = function(e){
    stop("The function u_to_param only supports families for which ktau_to_par is defined.
         This error likely is thrown, because simulate_non_simplified was called with
         the default argument for param_cond_func with a family,
         for which ktau_to_par is not implemented.")
  })
}

#' permutes indices of a rvine structure matrix so that on the antidiagonal
#' the values 1 to d (dimension of the matrix) are ordered from top right to bottom left
#' @param mat: The rvine structure matrix
permute_indices <- function(mat){
  antidiag_values <- rev(as_rvine_structure(mat)$order)
  n <- nrow(mat)
  permutation_dict <- list()
  for (i in 1:n){
    permutation_dict[[antidiag_values[i]]] <- i
  }
  struct_mat_new <- matrix(rep(0, length(mat)), nrow=nrow(mat))
  for (i in 1:n){
    for (j in 1:(n-i+1)){
      struct_mat_new[i,j] <- permutation_dict[[mat[i,j]]]
    }
  }
  return(list(struct_mat_new, permutation_dict, antidiag_values))
}

#' From the samples of the permuted structure, this function returns the samples from the original structure
#' @param simulations: Matrix with ncol(simulations) = length(original_order)
#' The samples from the reordered rvine structure (with antidiagonal values ordered
#' from bottom left to top right)
#' @param original_order: The elements on the antidiagonal from bottom left to top right in the original structure
permute_back <- function(simulations, original_order){
  result = simulations[,1:ncol(simulations)]
  for (i in 1:ncol(simulations)){
    result[,original_order[i]] <- simulations[,i]
  }
  return(result)
}

# function to compute the matrix \tilde{M}
get_max_matrix <- function(mat){
  num_cols <- ncol(mat)
  num_rows <- nrow(mat)
  max_mat <- matrix(rep(0, num_cols * num_rows), ncol=num_cols)
  for (j in 1:num_cols){
    max_in_col <- -1
    for (i in 1:(num_rows -j+1)){
      if(mat[i,j] > max_in_col){
        max_in_col <- mat[i,j]
      }
      max_mat[i,j] <- max_in_col
    }
  }
  return(max_mat)
}

#' Simulates a single sample,
#' This function is called by simulate_non_simplified, after modifying the input structure,
#' to be able to parallelize the procedure more easily for better performance
#' @param new_struct rvine structure, where the antidiagonal is ordered
#' increasingly from top right to bottom left.
#' @param M_mat the matrix \tilde{M}, with \tilde{M_{k,i}} = max_{j=1}^{k}{new_struct_{j,i}}
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
simulate_single_sample <- function(new_struct,
                                   M_mat,
                                   families,
                                   params,
                                   param_cond_funcs,
                                   rotations){
  d <- ncol(new_struct)
  w_row <- runif(d, min=0, max=1)
  u_row <- rep(0,d)
  # storage for the matrices V and V^2
  V <- matrix(rep(0,d^2), ncol=d)
  V_2 <- matrix(rep(0,d^2), ncol=d)
  V[1,d] <-w_row[d]
  # loop through the columns from right to left
  for (i in (d-1):1){
    V[d-i+1,i] <- w_row[i]
    # loop through the rows from the antidiagonal upwards
    for (k in (d-i):1){
      # In this iteration, we condition on m_{1,i},...,m_{k-1,i}
      # So for a non-simplified vine, the parameters can depend on
      # u_{m_{1,i}},..., u_{m_{k-1,i}}. These are already known at this point,
      # and they are stored in V_{1,{d-m_{1,i}+1}},...,V_{1,{d-m_{k-1,i}+1}}
      param <- params[[i]]
      if (k>1){
        cond_u_vals <- rep(0,k-1)
        for (temp in 1:(k-1)){
          cond_u_vals[temp] <- V[1,d-new_struct[temp,i]+1]
        }
        param <- param_cond_funcs[[k-1]][[i]](cond_u_vals, family=families[[k]][[i]])
      }
      cop <- bicop_dist(family=families[[k]][[i]],
                        rotation=rotations[[k]][[i]],
                        parameters = param)
      if(new_struct[k,i] == M_mat[k,i]){
        V[k,i] <- hbicop(c(V[k,d-M_mat[k,i]+1],
                           V[k+1,i]),
                         cop,
                         cond_var = 1, #i > m_{k,i} always holds for k <= d-i
                         inverse=TRUE)
      } else {
        V[k,i] <- hbicop(c(V_2[k, d-M_mat[k,i]+1],
                           V[k+1,i]),
                         cop,
                         cond_var = 1, #i > m_{k,i} always holds for k <= d-i
                         inverse=TRUE)
      }
      if(i>1){
        if(new_struct[k,i] == M_mat[k,i]){
          V_2[k+1,i] = hbicop(c(V[k,d-M_mat[k,i] + 1],
                                V[k,i]),
                              cop,
                              cond_var = 2, #i > m_{k,i} always holds for k <= d-i
                              inverse=FALSE)
        } else {
          V_2[k+1,i] = hbicop(c(V_2[k,d-M_mat[k,i] + 1],
                                V[k,i]),
                              cop,
                              cond_var=2, #i > m_{k,i} always holds for k <= d-i
                              inverse=FALSE)
        }
      }
    }
  }
  for(idx in 1:d){
    u_row[idx] <- V[1,d-idx+1]
  }
  return(u_row)
}

#' simulation of a non_simplified vine copula.
#' @param n_samples number of samples to simulate
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
simulate_non_simplified <- function(n_samples = 500,
                                    struct = matrix(c(1,1,1,
                                                      2,2,0,
                                                      3,0,0)
                                                    ,byrow=TRUE, ncol=3),
                                    families=list(list("gumbel", "gumbel"), list("gumbel")),
                                    params = list(c(2), c(1.3)),
                                    param_cond_funcs = list(list(u_to_param)),
                                    rotations = list(list(0,0),list(0))){
  # Reorder the indices, so they have the numbers 1 to d from top right to bottom left on the antidiagonal
  temp <- permute_indices(struct)
  new_struct <- temp[[1]]
  permutation <- temp[[2]]
  inverse_permutation <- temp[[3]]
  # Compute the matrix \tilde(M), with \tilde(M)_{k,i} = max(new_struct_{1,i},...,new_struct_{k,i})
  M_mat <- get_max_matrix(new_struct)
  # number of dimensions
  d <- ncol(new_struct)
  # matrix to save the simulations in
  u <- matrix(rep(0,d*n_samples),ncol=d)
  # loop once for every sample (For better performance in the future, maybe parallelize)
  for (sample in 1:n_samples){
    u_sample <- simulate_single_sample(new_struct,
                           M_mat,
                           families,
                           params,
                           param_cond_funcs,
                           rotations)
    u[sample,] <- u_sample
  }
  # permute the indices back to get a sample from the original rvine-structure
  result <- permute_back(u, inverse_permutation)
  return(result)
}

