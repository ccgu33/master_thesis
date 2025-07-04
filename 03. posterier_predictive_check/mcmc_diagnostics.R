library(rstan)
library(dplyr)

rds_base_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"

models <- list(
  list(
    name = "Complete Pooling",
    fit_file = file.path(rds_base_dir, "stan_fit_complete_pooling.rds")
  ),
  list(
    name = "No Pooling", 
    fit_file = file.path(rds_base_dir, "stan_fit_no_pooling.rds")
  ),
  list(
    name = "Partial Pooling Centered",
    fit_file = file.path(rds_base_dir, "stan_fit_partial_centered.rds")
  ),
  list(
    name = "Partial Pooling Non-Centered",
    fit_file = file.path(rds_base_dir, "stan_fit_partial_noncentered.rds")
  )
)


# Create results table
results <- data.frame(
  Model = character(),
  Parameters = integer(),
  Divergent_Transitions = integer(),
  Average_Rhat = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(models)) {
  model_info <- models[[i]]
  

  }
  
  fit <- readRDS(model_info$fit_file)
  
  # Get basic info
  summary_fit <- summary(fit)
  param_names <- rownames(summary_fit$summary)
  
  # Count parameters (exclude lp__)
  total_params <- sum(!grepl("lp__", param_names))
  
  # Get divergent transitions
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  # Get average Rhat (exclude lp__)
  rhat_values <- summary_fit$summary[, "Rhat"]
  rhat_clean <- rhat_values[!grepl("lp__", param_names)]
  rhat_clean <- rhat_clean[!is.na(rhat_clean) & is.finite(rhat_clean)]
  avg_rhat <- mean(rhat_clean)
  
  # Add to results
  results <- rbind(results, data.frame(
    Model = model_info$name,
    Parameters = total_params,
    Divergent_Transitions = divergent,
    Average_Rhat = round(avg_rhat, 4),
    stringsAsFactors = FALSE
  ))
}

# Display results
cat("MCMC DIAGNOSTICS SUMMARY:\n")
cat("=========================\n\n")

cat(sprintf("%-30s %12s %20s %15s\n", 
            "Model", "Parameters", "Divergent Transitions", "Average Rhat"))
cat(sprintf("%-30s %12s %20s %15s\n", 
            "-----", "----------", "-------------------", "------------"))

for(i in 1:nrow(results)) {
  row <- results[i,]
  cat(sprintf("%-30s %12d %20d %15.4f\n", 
              row$Model, row$Parameters, row$Divergent_Transitions, row$Average_Rhat))
}

cat("\n=== INTERPRETATION ===\n")
cat("Parameters: Total number of parameters in the model\n")
cat("Divergent Transitions: Number of divergent transitions (0 is ideal)\n")
cat("Average Rhat: Average R-hat across all parameters (1.00 is ideal, <1.01 is good)\n")

# Save results
write.csv(results, "simple_mcmc_diagnostics.csv", row.names = FALSE)
cat("\nResults saved to: simple_mcmc_diagnostics.csv\n")
