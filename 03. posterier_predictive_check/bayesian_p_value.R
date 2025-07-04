library(rstan)
library(dplyr)
library(tidyr)

rds_base_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"

#step 1: load the models
models <- list(
  list(
    name = "Complete Pooling",
    fit_file = file.path(rds_base_dir, "stan_fit_complete_pooling_normal.rds"),
    data_file = file.path(rds_base_dir, "stan_data_complete_pooling_normal.rds")
  ),
  list(
    name = "No Pooling", 
    fit_file = file.path(rds_base_dir, "stan_fit_no_pooling.rds"),
    data_file = file.path(rds_base_dir, "stan_data_no_pooling.rds")
  ),
  list(
    name = "Partial Pooling Centered",
    fit_file = file.path(rds_base_dir, "stan_fit_partial_centered.rds"),
    data_file = file.path(rds_base_dir, "stan_data_partial_centered.rds")
  ),
  list(
    name = "Partial Pooling Non-Centered",
    fit_file = file.path(rds_base_dir, "stan_fit_partial_noncentered.rds"),
    data_file = file.path(rds_base_dir, "stan_data_partial_noncentered.rds")
  )
)

extreme_low_threshold <- 0.05
extreme_high_threshold <- 0.95

#step 2: calculate the 7 p-values
results_list <- list()

for (i in 1:length(models)) {
  model_info <- models[[i]]
  fit <- readRDS(model_info$fit_file)
  stan_data <- readRDS(model_info$data_file)
  
  Y_obs_matrix <- stan_data$Y_obs # [10, 45] time x site matrix of observed coral cover
  Y_obs_present <- stan_data$Y_obs_present # [10, 45] binary matrix indicating which observations exist
  N_sites <- stan_data$N_sites #number of reef sites (45)
  N_years <- stan_data$N_years #number of time points (10)
  
  y_rep_array <- rstan::extract(fit, pars = "Y_rep")$Y_rep # [4000, 10, 45] iterations x time x site posterior predictive samples
  y_obs_vector <- as.vector(Y_obs_matrix) # [450] flattened vector of all observed values (complete data)
  n_observed <- length(y_obs_vector) # total number of data points (450)
  n_iter <- dim(y_rep_array)[1] # number of iterations (4000)
  y_rep_matrix_observed <- matrix(y_rep_array, nrow = n_iter, ncol = N_years * N_sites) # [4000, 450] all posterior samples reshaped
  
  y_rep_matrix_observed[y_rep_matrix_observed == -999] <- NA # replace missing value indicators with NA
  
    #mean
    # [450]observed data and then calculate the mean value, one single number
    T_obs_mean <- mean(y_obs_vector, na.rm = TRUE) 
    # each row is one posterior sample, and then calculate the mean value, [4000] vector
    T_rep_mean <- apply(y_rep_matrix_observed, 1, mean, na.rm = TRUE) 
    #how many out of 4000 posterior samples are greater than the observed mean?
    p_mean <- mean(T_rep_mean >= T_obs_mean)

    T_obs_sd <- sd(y_obs_vector, na.rm = TRUE)
    T_rep_sd <- apply(y_rep_matrix_observed, 1, sd, na.rm = TRUE)
    p_sd <- mean(T_rep_sd >= T_obs_sd)
    
    T_obs_min <- min(y_obs_vector, na.rm = TRUE)
    T_rep_min <- apply(y_rep_matrix_observed, 1, min, na.rm = TRUE)
    p_min <- mean(T_rep_min >= T_obs_min)
    
    T_obs_max <- max(y_obs_vector, na.rm = TRUE)
    T_rep_max <- apply(y_rep_matrix_observed, 1, max, na.rm = TRUE)
    p_max <- mean(T_rep_max >= T_obs_max)
    
    T_obs_q25 <- quantile(y_obs_vector, 0.25, na.rm = TRUE)
    T_rep_q25 <- apply(y_rep_matrix_observed, 1, quantile, 0.25, na.rm = TRUE)
    p_q25 <- mean(T_rep_q25 >= T_obs_q25)
    
    T_obs_q75 <- quantile(y_obs_vector, 0.75, na.rm = TRUE)
    T_rep_q75 <- apply(y_rep_matrix_observed, 1, quantile, 0.75, na.rm = TRUE)
    p_q75 <- mean(T_rep_q75 >= T_obs_q75)
    
    T_obs_median <- median(y_obs_vector, na.rm = TRUE)
    T_rep_median <- apply(y_rep_matrix_observed, 1, median, na.rm = TRUE)
    p_median <- mean(T_rep_median >= T_obs_median)
  
  p_values_calculated <- c(p_mean, p_sd, p_min, p_max, p_q25, p_q75, p_median) # [7] vector of all calculated p-values
  valid_p_values <- na.omit(p_values_calculated) # [7] vector with any NA values removed
  
        #percentage of extreme p-values
   pct_extreme <- mean(valid_p_values < extreme_low_threshold | valid_p_values > extreme_high_threshold) * 100 
  
   # Calculate average distance from 0.5 for all p-values
   avg_distance_from_05 <- mean(abs(valid_p_values - 0.5), na.rm = TRUE)
  
  results_list[[i]] <- data.frame(
    Model = model_info$name,
    p_mean,
    p_sd,
    p_min,
    p_max,
    p_q25,
    p_q75,
    p_median,
    pct_extreme,
    avg_distance_from_05,
    check.names = FALSE
  )
}

# step3: display the results
summary_df <- bind_rows(results_list) # data frame: [4 x 17] combined results from all models
print.data.frame(as.data.frame(summary_df), digits = 3, na.print = "NA", row.names = FALSE)


cat("SUMMARY TABLE: Average Distance from 0.5 and Extreme P-values\n")

final_table <- summary_df %>% # data frame: [4 x 4] 
  rowwise() %>%
  mutate(
    avg_pvalue = mean(c(p_mean, p_sd, p_min, p_max, 
                       p_q25, p_q75, p_median), na.rm = TRUE)
  ) %>%
  select(Model, avg_distance_from_05, pct_extreme, avg_pvalue) %>%
  arrange(avg_distance_from_05) %>%
  mutate(
    avg_distance_from_05 = round(avg_distance_from_05, 3),
    pct_extreme = round(pct_extreme, 3),
    avg_pvalue = round(avg_pvalue, 3)
  )

cat("\nFINAL TABLE:\n")
print(final_table, row.names = FALSE)
