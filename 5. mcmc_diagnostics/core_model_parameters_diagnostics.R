library(rstan)

# Load the model files
complete_pooling <- readRDS("/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds/stan_fit_complete_pooling.rds")
no_pooling <- readRDS("/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds/stan_fit_no_pooling.rds")
partial_centered <- readRDS("/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds/stan_fit_partial_centered.rds")
partial_noncentered <- readRDS("/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds/stan_fit_partial_noncentered.rds")

# Complete Pooling
summary1 <- summary(complete_pooling)
params1 <- rownames(summary1$summary)
core_params1 <- params1[!grepl("y_fit_tot|Y_rep|log_lik|lp__", params1)]
n_params1 <- length(core_params1)
divergent1 <- sum(get_divergent_iterations(complete_pooling))
rhat1 <- mean(summary1$summary[core_params1, "Rhat"], na.rm = TRUE)

# No Pooling
summary2 <- summary(no_pooling)
params2 <- rownames(summary2$summary)
core_params2 <- params2[!grepl("y_fit_tot|Y_rep|log_lik|lp__", params2)]
n_params2 <- length(core_params2)
divergent2 <- sum(get_divergent_iterations(no_pooling))
rhat2 <- mean(summary2$summary[core_params2, "Rhat"], na.rm = TRUE)

# Partial Pooling Centered
summary3 <- summary(partial_centered)
params3 <- rownames(summary3$summary)
core_params3 <- params3[!grepl("y_fit_tot|Y_rep|log_lik|lp__", params3)]
n_params3 <- length(core_params3)
divergent3 <- sum(get_divergent_iterations(partial_centered))
rhat3 <- mean(summary3$summary[core_params3, "Rhat"], na.rm = TRUE)

# Partial Pooling Non-Centered
summary4 <- summary(partial_noncentered)
params4 <- rownames(summary4$summary)
core_params4 <- params4[!grepl("y_fit_tot|Y_rep|log_lik|lp__", params4)]
n_params4 <- length(core_params4)
divergent4 <- sum(get_divergent_iterations(partial_noncentered))
rhat4 <- mean(summary4$summary[core_params4, "Rhat"], na.rm = TRUE)

# Print results
cat("Complete Pooling:", n_params1, "parameters,", divergent1, "divergent,", round(rhat1, 4), "rhat\n")
cat("No Pooling:", n_params2, "parameters,", divergent2, "divergent,", round(rhat2, 4), "rhat\n")
cat("Partial Pooling Centered:", n_params3, "parameters,", divergent3, "divergent,", round(rhat3, 4), "rhat\n")
cat("Partial Pooling Non-Centered:", n_params4, "parameters,", divergent4, "divergent,", round(rhat4, 4), "rhat\n")
