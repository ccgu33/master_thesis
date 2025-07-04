library(rstan)
library(loo)
library(dplyr)

rds_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/leave_one_out_validation"

#Step 1: load the models
model_files <- c(
  "stan_fit_complete_pooling.rds",
  "stan_fit_no_pooling.rds", 
  "stan_fit_partial_centered.rds",
  "stan_fit_partial_noncentered.rds"
)

# Model display names
model_names <- c(
  "Complete Pooling",
  "No Pooling",
  "Partial Pooling (Centered)",
  "Partial Pooling (Non-centered)"
)

#Step 2: load the data
data_file <- file.path(rds_dir, "stan_data_complete_pooling.rds")
stan_data <- readRDS(data_file)

# Get data dimensions
no_reps <- stan_data$N_sites
no_ts <- stan_data$N_years
n_data_points <- sum(stan_data$Y_obs_present)
t_gd <- stan_data$t_gd
ts <- stan_data$ts

# Storage for results
loo_results <- list()
log_lik_list <- list()

# Process each model
for (i in 1:length(model_files)) {
  
  fit <- readRDS(file.path(rds_dir, model_files[i]))
  
  # Extract samples
  y_fit_tot <- rstan::extract(fit, pars = "y_fit_tot")$y_fit_tot
  sigma_samples <- rstan::extract(fit, pars = "sigma")$sigma
  
  n_samples <- dim(y_fit_tot)[1]
  
  # Calculate log-likelihood matrix of data point d under sample s(full posterior distribution)
  log_lik <- matrix(0, nrow = n_samples, ncol = n_data_points) # rows: posterior samples, columns: data points, value: log-likelihood
  # For each posterior sample s
  for (s in 1:n_samples) {
    point_idx <- 1
    
    # For each site r
    for (r in 1:no_reps) {
      # For each time point t
      for (t in 1:no_ts) {
        if (stan_data$Y_obs_present[t, r] == 1) { #skip if the site is not present in the year(for the saprse data)
          # match prediction time to the closest time in the dense grid
          closest_t_idx <- which.min(abs(t_gd - ts[t]))
          
          # Get predicted and observed coral cover(actual parameter values)
          y_pred <- y_fit_tot[s, r, closest_t_idx] #sample s, site r, the matched time point
          y_obs <- stan_data$Y_obs[t, r] #observed coral cover
          
          # Handle different sigma structures across models
          if (length(dim(sigma_samples)) > 1) {
            # Site-specific sigma (no pooling, partial pooling models)
            sigma_r <- sigma_samples[s, r] 
          } else {
            # Shared sigma (complete pooling model)
            sigma_r <- sigma_samples[s]
          }
          
          # likelihood of data point d under the parameter values from posterior sample s
          log_lik[s, point_idx] <- dnorm(y_obs, mean = y_pred, sd = sigma_r, log = TRUE) 
          
          point_idx <- point_idx + 1 
        }
      }
    }
  }
  
  # Save log-likelihood matrix
  log_lik_list[[i]] <- log_lik
  
  # Compute LOO
  loo_results[[i]] <- loo(log_lik)
  print(loo_results[[i]])
  
  # Save individual LOO results
  saveRDS(loo_results[[i]], file.path(output_dir, paste0("loo_", gsub(" ", "_", tolower(model_names[i])), ".rds")))
}

# Save all log-likelihood matrices
saveRDS(log_lik_list, file.path(output_dir, "all_log_lik.rds"))

# Compare models
if (length(loo_results) > 1) {
  cat("\nModel comparison:\n")
  loo_comparison <- loo_compare(loo_results)
  print(loo_comparison)
  
  # Save comparison results
  saveRDS(loo_comparison, file.path(output_dir, "loo_comparison.rds"))
  
  # Create a summary table
  elpd_values <- numeric(length(loo_results))
  se_values <- numeric(length(loo_results))
  ploo_values <- numeric(length(loo_results))
  looic_values <- numeric(length(loo_results))
  
  # Extract values directly from LOO results matrix
  # [1] refers to the first model in loo_results list
  elpd_values[1] <- loo_results[[1]]$elpd_loo    # ELPD value
  se_values[1] <- loo_results[[1]]$se_elpd_loo   # ELPD standard error
  ploo_values[1] <- loo_results[[1]]$p_loo       # P_LOO value
  looic_values[1] <- loo_results[[1]]$looic      # LOOIC value
  
  elpd_values[2] <- loo_results[[2]]$elpd_loo    
  se_values[2] <- loo_results[[2]]$se_elpd_loo   
  ploo_values[2] <- loo_results[[2]]$p_loo       
  looic_values[2] <- loo_results[[2]]$looic      
  
  elpd_values[3] <- loo_results[[3]]$elpd_loo    
  se_values[3] <- loo_results[[3]]$se_elpd_loo   
  ploo_values[3] <- loo_results[[3]]$p_loo       
  looic_values[3] <- loo_results[[3]]$looic      
  
  elpd_values[4] <- loo_results[[4]]$elpd_loo    
  se_values[4] <- loo_results[[4]]$se_elpd_loo   
  ploo_values[4] <- loo_results[[4]]$p_loo       
  looic_values[4] <- loo_results[[4]]$looic      
  
  loo_summary <- data.frame(
    Model = model_names,
    ELPD = elpd_values,
    SE = se_values,
    P_LOO = ploo_values,
    LOOIC = looic_values
  )
  
  # Sort by ELPD (higher is better)
  loo_summary <- loo_summary %>% arrange(desc(ELPD))
  print(loo_summary)
  
  write.csv(loo_summary, file.path(output_dir, "loo_summary.csv"), row.names = FALSE)
  
 
 
}

