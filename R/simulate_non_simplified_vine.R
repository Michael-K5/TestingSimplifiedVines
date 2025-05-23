# This script contains functions to simulate from a Non-Simplified Vine Copula
library(rvinecopulib)

#' Computes the fisher Z transform (artanh)
#' @param tau: the value of Kendall's tau to transform to a real number.
#' Can also be a vector of taus.
#' @returns The fisher transformation(s) of tau. If tau is a number, this returns a number
#' if tau was a vector this returns a vector
fisher_z_transform = function(tau){
  return(atanh(tau))
}
#' Computes the inverse of the fisher Z transform, which is the tanh
#' @param T_val: The T-value that is supposed to be turned into a Kendall's tau again
#' @returns The inverse fisher transformation(s) of T_val. If T is a number, this returns a number
#' if T is a vector this returns a vector
inverse_fisher_transform <- function(T_val){
  # computes the inverse of the fisher Z transform (tanh)
  return(tanh(T_val))
}

#' Returns scaling_factor * tanh(T_val) + shift
#' @param T_val: The T-value that is supposed to be turned into a Kendall's tau again
#' @param scaling_factor: scale the result by this factor if desired. By default,
#' tanh returns values between -1 and 1. if a scaling_factor is specified, the result is
#' between - scaling factor and scaling factor
#' @param shift: shift the result if desired
#' @returns The inverse fisher transformation(s) of T_val. If T is a number, this returns a number
#' if T is a vector this returns a vector
scaled_tanh <- function(T_val, scaling_factor=1.0, shift=0.0){
  # computes the inverse of the fisher Z transform (tanh)
  return(scaling_factor * tanh(T_val) + shift)
}

#' maps a u value to the specific parameter for the third copula, by
#' finding a value on the fisher-z-transform-scale and mapping that back to a parameter
#' @param u: the samples u value (between 0 and 1)
#' @param family: the family of the third copula in the sample
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_param <- function(u, family="gaussian"){
  tryCatch({
    # tau should be in the range -0.92, +0.92 to avoid numerical problems later on
    tau_min <- -0.92
    tau_max <- 0.92
    T_min <- fisher_z_transform(tau_min)
    T_max <- fisher_z_transform(tau_max)
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
    stop(pasteo0("An error occurred: ",
    e,
    ". A typical cause of errors:
    The function u_to_param only supports families for which ktau_to_par is defined.
    This error likely is thrown, because simulate_non_simplified was called with
    the default argument for param_cond_funcs with a family,
    for which ktau_to_par is not implemented."))
  })
}

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function takes the dot product between a and u
#' and returns a linear function of that dot product.
#' @param a: a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
#' @param tau_lower: Lowest tau value, defaults to -0.92
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value, defaults to 0.92
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_param_linear <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      T_upper <- fisher_z_transform(tau_upper)
      T_lower <- fisher_z_transform(tau_lower)
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        T_lower <- fisher_z_transform(0.001)
      }
      T_val <- T_lower
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

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function takes the dot product between a and u
#' and returns a quadratic function of that dot product.
#' @param a: a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
#' @param tau_lower: Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_param_lin_quad <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      T_upper <- fisher_z_transform(tau_lower)
      T_lower <- fisher_z_transform(tau_upper)
      # Onepar families that cannot model negative dependence:
      # clayton, gumbel, joe
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        T_lower <- fisher_z_transform(0.001)
      }
      T_val <- T_lower
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

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function)
#' and then builds polynomial terms of degree 2, if dim(u) <=3, else a linear function of these
#' with a as parameters.
#' @param a: a vector of weights
#' @param tau_lower: Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_param_quadratic <- function(a, tau_lower=-0.92, tau_upper=0.92){
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
        # just apply linear function here.
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          The function u_to_param_non_lin only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function)
#' and then builds polynomial terms of degree 3, if dim(u) <=3, else a linear function of these
#' with a as parameters.
#' @param a: a vector of weights
#' @param tau_lower: Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper: Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @returns A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
u_to_param_cubic <- function(a, tau_lower=-0.92, tau_upper=0.92){
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
        # just apply linear function here.
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          The function u_to_param_non_lin only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' permutes indices of a rvine structure matrix so that on the antidiagonal
#' the values 1 to d (dimension of the matrix) are ordered from top right to bottom left
#' @param mat: The rvine structure matrix
#' @returns list with matrix as first entry, another list as second entry and a vector as third entry
#' First entry of result is the new structure matrix, second works as a dictionary, mapping the values
#' the last contains the values on the antidiagonal of the original matrix.
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
#' @returns result: the simulations, with the correct order.
permute_back <- function(simulations, original_order){
  result = simulations[,1:ncol(simulations)]
  for (i in 1:ncol(simulations)){
    result[,original_order[i]] <- simulations[,i]
  }
  return(result)
}

#' function to compute the matrix tilde(M)
#' @param mat: The input matrix
#' @returns the matrix, with entries tilde(m)_k,i = max_j=1,...,k (m_j,i)
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
#' @param M_mat the matrix tilde(M), with tilde(m)_k,i = max_{j=1}^{k}{new_struct_{j,i}}
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @returns u_row: a d-dimensional vector sampled from the non-simplified vine copula
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
#' For faster sampling please use simulate_non_simp_parallel instead.
#' @param n_samples number of samples to simulate
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @returns result: A nxd matrix, where each row contains one sample.
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
  # loop once for every sample
  prog_bar <- txtProgressBar(min=0,max=n_samples,title="Simulating...", style=3)
  for (sample in 1:n_samples){
    u_sample <- simulate_single_sample(new_struct,
                           M_mat,
                           families,
                           params,
                           param_cond_funcs,
                           rotations)
    u[sample,] <- u_sample
    setTxtProgressBar(prog_bar, value=sample, title="Simulating...")
  }
  # permute the indices back to get a sample from the original rvine-structure
  result <- permute_back(u, inverse_permutation)
  return(result)
}



#' parallel simulation of a non_simplified vine copula.
#' @param n_samples number of samples to simulate
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @returns result: A nxd matrix, where each row contains one sample.
simulate_non_simp_parallel <- function(n_samples = 500,
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

  # parallel computing
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # export the required functions for parallel processing
  parallel::clusterExport(cl,
                varlist = c("simulate_single_sample",
                            "u_to_param",
                            "bicop_dist",
                            "hbicop",
                            "ktau_to_par",
                            "fisher_z_transform",
                            "inverse_fisher_transform",
                            "scaled_tanh",
                            "struct",
                            "families",
                            "params",
                            "param_cond_funcs",
                            "rotations"),
                envir=environment())
  # pblapply to display a progress bar
  output_list <- pbapply::pblapply(1:n_samples, function(i) {
    simulate_single_sample(new_struct,
                           M_mat,
                           families,
                           params,
                           param_cond_funcs,
                           rotations)
  }, cl=cl)
  # stop the cluster, to free resources on computer
  parallel::stopCluster(cl)

  u <- do.call(rbind, output_list)
  # permute the indices back to get a sample from the original rvine-structure
  result <- permute_back(u, inverse_permutation)
  return(result)
}


# Plots
#' marginally normalized contour plots for copulas with different
#' families and one-dimensional conditioning values
#' @param cond_vals: The conditioning values
#' (plotted next to each other in columns)
#' @param family_vals: The names of the copula families in the rows
#' @param param_cond_func_1d: The parameter function to use
#' @param automatic_title: Whether to automatically create a title (TRUE or FALSE)
#' @param manual_title: Manually enter a title if desired. If this is not "",
#' then automatic_title is automatically set to FALSE
#' @param manual_subtitle: Manually enter a subtitle if desired. If this is "",
#' then an automatic subtitle is set,
#' if either automatic_title is TRUE or manual_title is not "".
plot_contours_1d <- function(
    cond_vals=c(0.2,0.4,0.6,0.8),
    family_vals=c("frank","gaussian","gumbel", "clayton"),
    param_cond_func_1d = u_to_param_linear(c(1)),
    automatic_title=TRUE,
    manual_title="",
    manual_subtitle=""){
  n <- length(cond_vals)
  m <- length(family_vals)
  # leave enough space for a title on top of the plots, if desired
  space_on_top <- 1.5
  if(automatic_title || manual_title != ""){
    space_on_top <- 4
  }
  par(mfrow=c(m, n), mar=c(1,1,1,1), oma=c(1.5,1.5,space_on_top,1.5))
  for(i in  1:m){
    for(j in 1:n){
      temp_param <- param_cond_func_1d(cond_vals[j], family=family_vals[i])
      bicop_cont_dist <- bicop_dist(
        family=family_vals[i],
        rotation=0,
        parameters = temp_param)
      contour(bicop_cont_dist,margins="norm",
              main=paste(family_vals[i], round(cond_vals[j],2)),
              drawlabels=FALSE,
              axes=FALSE)
      box()
    }
  }
  if(automatic_title && manual_title == ""){
    mtext("Marginally normalized contour plots", outer = TRUE, cex = 1, line = 2, font=2)
    if(manual_subtitle==""){
      mtext("Different families and conditioning values, see titles", outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle, outer = TRUE, cex = 0.8, line = 1)
    }
  }
  if(manual_title != ""){
    mtext(manual_title, outer=TRUE, cex=1, line=2, font=2)
    if(manual_subtitle==""){
      mtext("Different families and conditioning values, see titles", outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle, outer = TRUE, cex = 0.8, line = 1)
    }
  }
}

#' marginally normalized contour plots for a copula with family family_name
#' with different two-dimensional conditioning values
#' @param cond_vals_1: The conditioning values of the first variable
#' (plotted next to each other in columns)
#' @param cond_vals_2: The conditioning values of the second variable
#' (plotted in the rows)
#' @param family_name: The name of the copula family to use
#' @param param_cond_func: The parameter function to use
#' @param automatic_title: Whether an automatically created title
#' should be plotted.Needs to be TRUE or FALSE.
#' @param manual_title: Manual Title, if desired. If this is not "", then
#' automatic_title is automatically set to FALSE.
#' @param manual_subtitle: Manual Subtitle, if desired. If this is "",
#' then an automatic subtitle is set,
#' if either automatic_title is TRUE or manual_title is not "".
plot_contours_2d <- function(
    cond_vals_1 = c(0.2,0.4,0.6,0.8),
    cond_vals_2 = c(0.2,0.4,0.6,0.8),
    family_name="frank",
    param_cond_func_2d = u_to_param_linear(c(0.7,0.3)),
    automatic_title = TRUE,
    manual_title="",
    manual_subtitle=""){
  m <- length(cond_vals_1)
  n <- length(cond_vals_2)
  # leave enough space for a title on top of the plots, if desired
  space_on_top <- 1.5
  if(automatic_title || manual_title != ""){
    space_on_top <- 4
  }
  par(mfrow=c(m, n), mar=c(1,1,1,1), oma=c(1.5,1.5,space_on_top,1.5))
  for(i in  1:m){
    for(j in 1:n){
      temp_param <- param_cond_func_2d(c(cond_vals_1[i], cond_vals_2[j]), family=family_name)
      bicop_cont_dist <- bicop_dist(
        family=family_name,
        rotation=0,
        parameters = temp_param)
      contour(bicop_cont_dist,margins="norm",
              main=paste("1:", cond_vals_1[i], "2:", cond_vals_2[j]),
              drawlabels=FALSE,
              axes=FALSE)
      box()
    }
  }
  if(automatic_title && manual_title == ""){
    mtext(paste0("Marginally normalized contour plots for the ", family_name, " copula"),
          outer = TRUE, cex = 1, line = 2, font=2)
    if(manual_subtitle==""){
      mtext("Two conditioning values, see plot titles",
            outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle,
            outer = TRUE, cex = 0.8, line = 1)
    }
  }
  if(manual_title != ""){
    mtext(manual_title, outer=TRUE, cex=1, line=2, font=2)
    if(manual_subtitle==""){
      mtext("Two conditioning values, see plot titles",
            outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle,
            outer = TRUE, cex = 0.8, line = 1)
    }
  }
}



#' Function with similar functionality as the pairs_copula_data function of
#' the rvinecopulib package. Plots histograms of df on the diagonal,
#' scatter plots on the upper triangle and
#' marginally normalized contours on the lower triangle.
#' @param u_data Data to plot, each dimension needs to be between 0 and 1
#' @param max_samples Maximum number of samples to consider for the plot
#' If -1, then all samples are considered.
#' If >0, then a random sample of max_samples is drawn and considered for the plots
#' @param plot_emp_ktau Whether to plot the empirical kendall's tau (TRUE or FALSE)
pairs_copula_data_custom <- function(u_data, max_samples = -1, plot_emp_ktau=FALSE) {
  df <- u_data
  if (is.matrix(df)) {
    df <- as.data.frame(df)
  }
  if(!is.data.frame(df) && !is.matrix(df)){
    stop("u_data needs to be either a dataframe or a matrix")
  }
  n <- ncol(df)
  if(max_samples > 0 && max_samples < n){
    to_plot_idx <- sample(n,floor(max_samples))
    df <- u_data[to_plot_idx,]
  }
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
        if(plot_emp_ktau){
          emp_k_tau <- cor(df[[j]], df[[i]], method="kendall")
          mtext(round(emp_k_tau,2), line=-3, col="blue",font=2)
        }
        box()
      } else {
        # Contour plot using 2D density estimation
        x <- df[[j]]
        y <- df[[i]]
        # convert to R using normal quantile function (yields easier to read contours)
        x <- qnorm(x)
        y <- qnorm(y)
        kde <- MASS::kde2d(x, y, n = 30, lims = c(-3,3,-3,3)) #lims = c(0, 1, 0, 1)
        contour(kde, nlevels=8,drawlabels = FALSE, xlab = "", ylab = "", axes = FALSE)
        box()
      }
    }
  }
}

