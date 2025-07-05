library(rstan)
library(shinystan)
rds_dir <- '/Users/chengu/Documents/master thesis/Bayesian-inference-of-coral-bleaching-dynamics/01. pooling_methods/rds'

launch_shinystan_for_model <- function(model_name) {
  fit_file <- file.path(rds_dir, paste0("stan_fit_", model_name, ".rds"))
  fit <- readRDS(fit_file)
  sso <- as.shinystan(fit)
  launch_shinystan(sso, launch.browser = FALSE)
}

launch_shinystan_for_model("partial_noncentered")