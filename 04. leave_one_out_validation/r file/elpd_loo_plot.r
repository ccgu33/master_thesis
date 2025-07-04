library(ggplot2)
library(dplyr)
library(bayesplot)
library(gridExtra)
library(loo)
library(patchwork) # For combining plots

rds_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
loo_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/leave_one_out_validation"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/plots/leave_one_out_validation"



# Model names
model_names <- c(
  "Complete Pooling",
  "No Pooling",
  "Partial Pooling (Centered)",
  "Partial Pooling (Non-centered)"
)

# Load LOO results
loo_results <- list()
for (i in 1:length(model_names)) {
  loo_file <- file.path(loo_dir, paste0("loo_", gsub(" ", "_", tolower(model_names[i])), ".rds"))
  loo_results[[i]] <- readRDS(loo_file)
}

# 1. Combined Pareto k diagnostic plot for all models
# Create a data frame with all Pareto k values
all_k_data <- data.frame()

for (i in 1:length(loo_results)) {
  k_values <- loo_results[[i]]$diagnostics$pareto_k
  
  # Add to combined data frame
  model_k_data <- data.frame(
    index = 1:length(k_values),
    k = k_values,
    model = model_names[i]
  )
  
  all_k_data <- rbind(all_k_data, model_k_data) #index, elpd, model name(combined all the models)
}

# Create faceted plot
p_k_combined <- ggplot(all_k_data, aes(x = index, y = k, color = model)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0.7, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "orange", linetype = "dashed") +
  facet_wrap(~ model, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Pareto k Diagnostic for All Models",
    x = "Data Point Index",
    y = "Pareto k"
  ) +
  theme_minimal(base_size = 12) + #get help from chatgpt
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

# Save the combined Pareto k plot
ggsave(file.path(output_dir, "combined_pareto_k_plot.png"), p_k_combined, width = 10, height = 8)

# 2. Combined pointwise ELPD plot for all models
all_elpd_data <- data.frame()

for (i in 1:length(loo_results)) {
  # Get pointwise elpd values
  pointwise_elpd <- loo_results[[i]]$pointwise[,"elpd_loo"]
  
  # Add to combined data frame
  model_elpd_data <- data.frame(
    index = 1:length(pointwise_elpd),
    elpd = pointwise_elpd,
    model = model_names[i]
  )
  
  all_elpd_data <- rbind(all_elpd_data, model_elpd_data)
}

# Create faceted plot
p_elpd_combined <- ggplot(all_elpd_data, aes(x = index, y = elpd, color = model)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~ model, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Pointwise ELPD for All Models",
    subtitle = "Higher values indicate better fit for individual observations",
    x = "Data Point Index",
    y = "ELPD"
  ) +
  theme_minimal(base_size = 12) + #get help from chatgpt
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

# Save the combined ELPD plot
ggsave(file.path(output_dir, "combined_pointwise_elpd_plot.png"), p_elpd_combined, width = 10, height = 8)

# 3. Create a single combined diagnostic plot with both visualizations (get help from chatgpt)
p_combined <- p_k_combined / p_elpd_combined
ggsave(file.path(output_dir, "combined_diagnostics.png"), p_combined, width = 12, height = 14)

cat("Combined diagnostic plots created and saved to:", output_dir, "\n")
