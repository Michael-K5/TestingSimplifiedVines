# Evaluation of the experiments
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
# Create a vector of file names
all_filenames <- paste0("results/20250617OneParNSamples10000Repl", 1:10, ".csv")
print(all_filenames)
# sum(is.na(all_results))
#quantile(all_results["MC_Integral"], probs=c(0.05,0.1,0.9,0.95), na.rm=TRUE)

# Read and combine all CSVs into one data frame
all_results <- bind_rows(lapply(all_filenames, read_csv))
all_results <- all_results %>%
  mutate(param_cond_func_idx = recode(param_cond_func_idx,
                                      `1` = "linear",
                                      `2` = "quadratic",
                                      `3` = "cubic")) %>%
  rename(param_cond_func = param_cond_func_idx) %>%
  rename(tau_max = tau_max)
# all_results
#all_results[1:5,"lin-q0.05geq0"]
summary_df <- all_results %>%
  group_by(dim, tau_max, param_cond_func) %>%
  summarise(
    onepar0.05_avg = mean(`onepar-q0.05geq0`, na.rm = TRUE),
    onepar0.05_std   = sd(`onepar-q0.05geq0`, na.rm = TRUE),
    onepar0.1_avg = mean(`onepar-q0.1geq0`, na.rm = TRUE),
    onepar0.1_std   = sd(`onepar-q0.1geq0`, na.rm = TRUE),
    onepar0.95_avg = mean(`onepar-q0.95leq0`, na.rm = TRUE),
    onepar0.95_std   = sd(`onepar-q0.95leq0`, na.rm = TRUE),
    onepar0.9_avg = mean(`onepar-q0.9leq0`, na.rm = TRUE),
    onepar0.9_std   = sd(`onepar-q0.9leq0`, na.rm = TRUE),
    nonpar0.05_avg = mean(`nonpar-q0.05geq0`, na.rm = TRUE),
    nonpar0.05_std   = sd(`nonpar-q0.05geq0`, na.rm = TRUE),
    nonpar0.1_avg = mean(`nonpar-q0.1geq0`, na.rm = TRUE),
    nonpar0.1_std   = sd(`nonpar-q0.1geq0`, na.rm = TRUE),
    nonpar0.95_avg = mean(`nonpar-q0.95leq0`, na.rm = TRUE),
    nonpar0.95_std   = sd(`nonpar-q0.95leq0`, na.rm = TRUE),
    nonpar0.9_avg = mean(`nonpar-q0.9leq0`, na.rm = TRUE),
    nonpar0.9_std   = sd(`nonpar-q0.9leq0`, na.rm = TRUE),
    lin0.05_avg = mean(`lin-q0.05geq0`, na.rm = TRUE),
    lin0.05_std   = sd(`lin-q0.05geq0`, na.rm = TRUE),
    lin0.1_avg = mean(`lin-q0.1geq0`, na.rm = TRUE),
    lin0.1_std   = sd(`lin-q0.1geq0`, na.rm = TRUE),
    lin0.95_avg = mean(`lin-q0.95leq0`, na.rm = TRUE),
    lin0.95_std   = sd(`lin-q0.95leq0`, na.rm = TRUE),
    lin0.9_avg = mean(`lin-q0.9leq0`, na.rm = TRUE),
    lin0.9_std   = sd(`lin-q0.9leq0`, na.rm = TRUE),
    .groups = 'drop'  # ungroup after summarising
  )

summary_df_05 <- all_results %>%
  group_by(dim, tau_max, param_cond_func) %>%
  summarise(
    onepar0.05_avg = round(mean(`onepar-q0.05geq0`, na.rm = TRUE),3),
    onepar0.05_std   = round(sd(`onepar-q0.05geq0`, na.rm = TRUE),3),
    nonpar0.05_avg = round(mean(`nonpar-q0.05geq0`, na.rm = TRUE),3),
    nonpar0.05_std   = round(sd(`nonpar-q0.05geq0`, na.rm = TRUE),3),
    lin0.05_avg = round(mean(`lin-q0.05geq0`, na.rm = TRUE),3),
    lin0.05_std   = round(sd(`lin-q0.05geq0`, na.rm = TRUE),3),
    .groups = 'drop'  # ungroup after summarising
  )
summary_df_05
knitr::kable(summary_df_05, format="latex", booktabs="TRUE")

summary_df_10 <- all_results %>%
  group_by(dim, tau_max, param_cond_func) %>%
  summarise(
    onepar0.10_avg = round(mean(`onepar-q0.1geq0`, na.rm = TRUE),3),
    onepar0.10_std   = round(sd(`onepar-q0.1geq0`, na.rm = TRUE),3),
    nonpar0.10_avg = round(mean(`nonpar-q0.1geq0`, na.rm = TRUE),3),
    nonpar0.10_std   = round(sd(`nonpar-q0.1geq0`, na.rm = TRUE),3),
    lin0.10_avg = round(mean(`lin-q0.1geq0`, na.rm = TRUE),3),
    lin0.10_std   = round(sd(`lin-q0.1geq0`, na.rm = TRUE),3),
    .groups = 'drop'  # ungroup after summarising
  )
summary_df_10
knitr::kable(summary_df_10, format="latex", booktabs="TRUE")

summary_df_90 <- all_results %>%
  group_by(dim, tau_max, param_cond_func) %>%
  summarise(
    onepar0.90_avg = round(mean(`onepar-q0.9leq0`, na.rm = TRUE),3),
    onepar0.90_std   = round(sd(`onepar-q0.9leq0`, na.rm = TRUE),3),
    nonpar0.90_avg = round(mean(`nonpar-q0.9leq0`, na.rm = TRUE),3),
    nonpar0.90_std   = round(sd(`nonpar-q0.9leq0`, na.rm = TRUE),3),
    lin0.90_avg = round(mean(`lin-q0.9leq0`, na.rm = TRUE),3),
    lin0.90_std   = round(sd(`lin-q0.9leq0`, na.rm = TRUE),3),
    .groups = 'drop'  # ungroup after summarising
  )
summary_df_90
knitr::kable(summary_df_90, format="latex", booktabs="TRUE")

summary_df_95 <- all_results %>%
  group_by(dim, tau_max, param_cond_func) %>%
  summarise(
    onepar0.90_avg = round(mean(`onepar-q0.95leq0`, na.rm = TRUE),3),
    onepar0.90_std   = round(sd(`onepar-q0.95leq0`, na.rm = TRUE),3),
    nonpar0.90_avg = round(mean(`nonpar-q0.95leq0`, na.rm = TRUE),3),
    nonpar0.90_std   = round(sd(`nonpar-q0.95leq0`, na.rm = TRUE),3),
    lin0.90_avg = round(mean(`lin-q0.95leq0`, na.rm = TRUE),3),
    lin0.90_std   = round(sd(`lin-q0.95leq0`, na.rm = TRUE),3),
    .groups = 'drop'  # ungroup after summarising
  )
summary_df_95
knitr::kable(summary_df_95, format="latex", booktabs="TRUE")


prep_data_for_plot <- function(df, q_reg_type = "onepar"){
  # Step 1: Pivot the data longer
  df_long <- df %>%
    pivot_longer(
      cols = starts_with(paste0(q_reg_type, "-q")),
      names_to = "quantile",
      values_to = "value"
    ) %>%
    # Step 2: Extract numeric quantile from the column name
    mutate(
      quantile = as.numeric(sub(paste0(q_reg_type, "-q([0-9.]+)geq0"), "\\1", quantile))
    ) %>%
    filter(quantile <= 0.3)

  df_long <- df_long %>%
    group_by(dim, tau_max, param_cond_func, quantile) %>%
    mutate(sim_id = row_number()) %>%
    ungroup()
  return(df_long)
}
df_long <- prep_data_for_plot(df= all_results, q_reg_type="onepar")

df_long_dim3 <- df_long %>%
  filter(dim == 3)

df_long_dim5 <- df_long %>%
  filter(dim ==5)

plot_colored_dashed_lines <- function(df){
  p <- ggplot(df, aes(x = quantile, y = value,
                           color = factor(tau_max),
                           linetype = factor(param_cond_func),
                           group = interaction(dim, tau_max, param_cond_func, sim_id))) +
    geom_line(alpha = 0.6) +  # alpha < 1 to make overlapping lines visible
    facet_wrap(~ dim, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("0.3" = "black", "0.6" = "blue", "0.9" = "orange")) +
    scale_linetype_manual(values = c("linear" = "solid", "quadratic" = "dashed", "cubic" = "dotted")) +
    labs(
      title = "Test Results",
      x = "Quantile level",
      y = "Number of Points",
      color = "Tau Max",
      linetype = "Param Cond Func"
    ) +
    theme_minimal(base_size = 14)
  return(p)
}
p <- plot_colored_dashed_lines(df_long_dim3)
p
ggsave(
  filename = paste0("plot_dim_", 3, ".png"),
  plot = p,
  width = 8,
  height = 4,
  dpi = 300,
  bg="white"
)
plot_colored_dashed_lines(df_long_dim5)


# Plot mean and confidence bands
plot_mean_conf <- function(df_long){
  # Calculate summary stats: mean and standard error
  summary_df <- df_long %>%
    group_by(dim, tau_max, param_cond_func, quantile) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  p <- ggplot() +
    # Individual simulation lines (light and transparent)
    geom_line(data = df_long,
              aes(x = quantile, y = value,
                  color = factor(tau_max),
                  linetype = factor(param_cond_func),
                  group = interaction(dim, tau_max, param_cond_func, sim_id)),
              alpha = 0.3) +
    # Mean lines
    geom_line(data = summary_df,
              aes(x = quantile, y = mean_value,
                  color = factor(tau_max),
                  linetype = factor(param_cond_func)),
              linewidth = 1) +
    # Confidence ribbon
    geom_ribbon(data = summary_df,
                aes(x = quantile,
                    ymin = mean_value - 1.96 * se_value,
                    ymax = mean_value + 1.96 * se_value,
                    fill = factor(tau_max),
                    group = interaction(tau_max, param_cond_func)),
                alpha = 0.2,
                color = NA) +
    facet_wrap(~ dim, ncol = 1, scales = "free_y", labeller= label_both) +
    scale_color_manual(values = c("0.3" = "black", "0.6" = "blue", "0.9" = "orange")) +
    scale_fill_manual(values = c("0.3" = "black", "0.6" = "blue", "0.9" = "orange")) +
    scale_linetype_manual(values = c("linear" = "solid", "quadratic" = "dashed", "cubic" = "dotted")) +
    labs(
      x = "Quantile level",
      y = "Number of Points",
      color = "Tau Max",
      fill = "Tau Max",
      linetype = "Param Cond Func"
    ) +
    theme_minimal(base_size = 14)
  return(p)
}

p3 <- plot_mean_conf(df_long_dim3)
p3
p5 <- plot_mean_conf(df_long_dim5)
p5



prep_data_for_plot_upper <- function(df, q_reg_type = "onepar"){
  # Step 1: Pivot the data longer
  df_long <- df %>%
    pivot_longer(
      cols = starts_with(paste0(q_reg_type, "-q")),
      names_to = "quantile",
      values_to = "value"
    ) %>%
    # Step 2: Extract numeric quantile from the column name
    mutate(
      quantile = as.numeric(sub(paste0(q_reg_type, "-q([0-9.]+)leq0"), "\\1", quantile))
    ) %>%
    filter(quantile >= 0.7)

  df_long <- df_long %>%
    group_by(dim, tau_max, param_cond_func, quantile) %>%
    mutate(sim_id = row_number()) %>%
    ungroup()
  return(df_long)
}
df_long_upper <- prep_data_for_plot_upper(df= all_results, q_reg_type="nonpar")
df_long_upper_dim3 <- df_long_upper %>%
  filter(dim==3)
df_long_upper_dim5 <- df_long_upper %>%
  filter(dim==5)
plot_colored_dashed_lines(df_long_upper_dim3)
plot_colored_dashed_lines(df_long_upper_dim5)
plot_mean_conf(df_long_upper_dim3)
plot_mean_conf(df_long_upper_dim5)


# All relevant:
library(patchwork)
q_reg_types = c("onepar", "nonpar", "lin")
for(q_reg_type in q_reg_types){
  df_long <- prep_data_for_plot(df= all_results, q_reg_type=q_reg_type)
  df_long_dim3 <- df_long %>%
    filter(dim == 3)
  df_long_dim5 <- df_long %>%
    filter(dim ==5)
  p3_lower <- plot_mean_conf(df_long_dim3)
  p5_lower <- plot_mean_conf(df_long_dim5)
  df_long_upper <- prep_data_for_plot_upper(df= all_results, q_reg_type=q_reg_type)
  df_long_upper_dim3 <- df_long_upper %>%
    filter(dim==3)
  df_long_upper_dim5 <- df_long_upper %>%
    filter(dim==5)
  p3_upper <- plot_mean_conf(df_long_upper_dim3)
  p5_upper <- plot_mean_conf(df_long_upper_dim5)
  q_reg_print_name <- "One Parameter D-Vine"
  if(q_reg_type == "nonpar"){
    q_reg_print_name <- "Nonparametric D-Vine"
  } else if(q_reg_type =="lin"){
    q_reg_print_name <- "Linear"
  }
  joint_plot <- (p3_lower + p3_upper) / (p5_lower + p5_upper) +
    plot_annotation(
      title = paste0("Test Results for ",q_reg_print_name," Quantile Regression"),
    )
  ggsave(
    filename = paste0(q_reg_type, "_plot.png"),
    plot = joint_plot,
    width = 12,
    height = 8,
    dpi = 300,
    bg="white"
  )
}





#TEMP
# library(stringr)
# fix_df_names <- function(i){
#   filename <- paste0("results/20250617OneParnSamples10000Repl",i,".csv")
#   temp_df <- read_csv(filename)
#   temp_df <- temp_df %>%
#     rename_with(
#       .cols = matches("^onepar-q[0-9.]+geq0$"),
#       .fn = function(colnames) {
#         sapply(colnames, function(colname) {
#           qval <- as.numeric(str_match(colname, "^onepar-q([0-9.]+)geq0$")[, 2])
#           if (!is.na(qval) && qval > 0.7) {
#             str_replace(colname, "geq", "leq")
#           } else {
#             colname
#           }
#         })
#       }
#     )
#   temp_df <- temp_df %>%
#     rename_with(
#       .cols = matches("^nonpar-q[0-9.]+geq0$"),
#       .fn = function(colnames) {
#         sapply(colnames, function(colname) {
#           qval <- as.numeric(str_match(colname, "^nonpar-q([0-9.]+)geq0$")[, 2])
#           if (!is.na(qval) && qval > 0.7) {
#             str_replace(colname, "geq", "leq")
#           } else {
#             colname
#           }
#         })
#       }
#     )
#   temp_df <- temp_df %>%
#     rename_with(
#       .cols = matches("^lin-q[0-9.]+geq0$"),
#       .fn = function(colnames) {
#         sapply(colnames, function(colname) {
#           qval <- as.numeric(str_match(colname, "^lin-q([0-9.]+)geq0$")[, 2])
#           if (!is.na(qval) && qval > 0.7) {
#             str_replace(colname, "geq", "leq")
#           } else {
#             colname
#           }
#         })
#       }
#     )
#   result_filename <- paste0("results/20250617OneParNSamples10000Repl",i,"fixed.csv")
#   write_csv(temp_df, result_filename)
# }
# for(i in 1:10){
#   fix_df_names(i)
# }
#
# compare_dataframes_fuzzy <- function(df1, df2, tol = 1e-3) {
#     if (!all(dim(df1) == dim(df2))) {
#       return(FALSE)
#     }
#
#     mat1 <- as.matrix(df1)
#     mat2 <- as.matrix(df2)
#
#     for (i in seq_len(nrow(mat1))) {
#       for (j in seq_len(ncol(mat1))) {
#         val1 <- mat1[i, j]
#         val2 <- mat2[i, j]
#
#         # Both NA → equal
#         if (is.na(val1) && is.na(val2)) next
#
#         # Only one NA → not equal
#         if (is.na(val1) || is.na(val2)) return(FALSE)
#
#         # Numeric comparison with tolerance
#         if (is.numeric(val1) && is.numeric(val2)) {
#           if (abs(val1 - val2) > tol) return(FALSE)
#         } else {
#           # Non-numeric comparison must be exact
#           if (val1 != val2) return(FALSE)
#         }
#       }
#     }
#
#     return(TRUE)
# }
# vals <- list()
# for(i in 1:10){
#   filename_old <- paste0("results/20250617OneParnSamples10000Repl",i,".csv")
#   filename_fixed <- paste0("results/20250617OneParNSamples10000Repl",i,"fixed.csv")
#   df_old <- read_csv(filename_old)
#   df_fixed <- read_csv(filename_fixed)
#   vals[[i]] <- compare_dataframes_fuzzy(df_old, df_fixed)
# }
# print(vals)

