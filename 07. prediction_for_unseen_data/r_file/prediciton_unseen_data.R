FORGETTING_FACTOR <- 0.9

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(purrr)

rds_base_dir <- '/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds'
raw_data_file <- '/Users/chengu/Documents/master thesis/new_reef_data.csv'
fit_file <- file.path(rds_base_dir, "stan_fit_no_pooling.rds")
data_file <- file.path(rds_base_dir, "stan_data_no_pooling.rds")
plot_output_dir <- file.path('/Users/chengu/Documents/master thesis/final_final_final_code/plots/prediction_for_unseen_data', "validation_plots_from_fit_model")

site_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
START_YEAR <- 1995
END_YEAR_TRAIN <- 2004
END_YEAR_VALID <- 2007

fit <- readRDS(fit_file)
stan_data <- readRDS(data_file)

reef_names <- readRDS(file.path(rds_base_dir, "reef_id_mapping.rds"))$REEF_NAME

site_map <- data.frame(
  SiteID = 1:45,
  REEF_NAME = rep(reef_names, each = 3),
  SiteWithinReef = rep(1:3, 15),
  ReefID = rep(1:15, each = 3)
)

raw_data <- read.csv(raw_data_file, sep = ';', header = TRUE, stringsAsFactors = FALSE) %>%
  select(REEF_NAME, SITE_NO, SAMPLE_DATE, GROUP_CODE, COVER)

observed_data_long <- raw_data %>%
  filter(GROUP_CODE == 'Hard Coral') %>%
  mutate(
    YEAR = as.numeric(substr(SAMPLE_DATE, 1, 4)),
    Cover = as.numeric(gsub(',', '.', COVER)) / 100
  ) %>%
  filter(YEAR >= START_YEAR, YEAR <= END_YEAR_VALID) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarize(Cover = mean(Cover, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(site_map, by = c("REEF_NAME" = "REEF_NAME", "SITE_NO" = "SiteWithinReef")) %>%
  mutate(Period = ifelse(YEAR <= END_YEAR_TRAIN, "Training", "Validation"))

# Extract parameters
params <- rstan::extract(fit, pars = c("theta", "y0_init_C"))

# Generate predictions
prediction_times <- seq(0, END_YEAR_VALID - START_YEAR, by = 0.1)
set.seed(123)
n_samples <- min(1000, dim(params$theta)[1])
sample_indices <- sample(n_samples)

# Solve ODEs
all_trajectories <- list()

for (i in seq_along(sample_indices)) {
  sample_idx <- sample_indices[i]
  
  for (s in 1:stan_data$N_sites) {
    theta_s <- pmax(1e-6, params$theta[sample_idx, , s])
    N0 <- params$y0_init_C[sample_idx, s]
    y0 <- c(C = pmax(1e-6, pmin(1-1e-6, N0^2)), 
            B = pmax(1e-6, pmin(1-1e-6, N0 * (1-N0))))
    
    # ODE function inline
    coral_ode <- function(t, y, parms) {
      C <- pmax(0, pmin(1, y[1]))
      B <- pmax(0, pmin(1, y[2]))
      available <- pmax(0, 1 - C - B)
      dCdt <- parms[1] * C * available - parms[2] * C + parms[3] * B
      dBdt <- parms[2] * C - parms[3] * B - parms[4] * B
      list(c(dCdt, dBdt))
    }
    
    ode_result <- deSolve::ode(y = y0, times = prediction_times, func = coral_ode, parms = theta_s)
    
    all_trajectories[[length(all_trajectories) + 1]] <- data.frame(
      Sample = sample_idx, 
      SiteID = s, 
      Time = ode_result[,1],
      TotalCover = pmin(1, pmax(0, rowSums(ode_result[,2:3])))
    )
  }
}

# Summarize predictions
predicted_summary_long <- bind_rows(all_trajectories) %>%
  group_by(SiteID, Time) %>%
  summarize(
    LowerCI = quantile(TotalCover, 0.025, na.rm = TRUE),
    MedianCover = quantile(TotalCover, 0.5, na.rm = TRUE),
    UpperCI = quantile(TotalCover, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(YEAR = Time + START_YEAR) %>% 
  left_join(site_map, by = "SiteID")

# Set factor levels for both datasets at once
predicted_summary_long$REEF_NAME <- factor(predicted_summary_long$REEF_NAME, levels = reef_names)
observed_data_long$REEF_NAME <- factor(observed_data_long$REEF_NAME, levels = reef_names)

# Create validation plot
validation_plot <- ggplot() +
  geom_ribbon(data = predicted_summary_long, 
              aes(x = YEAR, ymin = LowerCI, ymax = UpperCI, 
                  group = factor(SiteWithinReef), fill = factor(SiteWithinReef)), 
              alpha = 0.3) +
  geom_line(data = predicted_summary_long, 
            aes(x = YEAR, y = MedianCover, 
                group = factor(SiteWithinReef), color = factor(SiteWithinReef)), 
            linewidth = 1) +
  geom_point(data = observed_data_long, 
             aes(x = YEAR, y = Cover, color = factor(SITE_NO), shape = Period), 
             size = 2.5, alpha = 0.8) +
  geom_vline(xintercept = END_YEAR_TRAIN + 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_shape_manual(name = "Period", values = c("Training" = 16, "Validation" = 17)) +
  facet_wrap(~ REEF_NAME, scales = "free_y", ncol = 3) +
  labs(title = "Predicted vs. Observed Coral Cover (1995-2007)",
       subtitle = "Lines are median predictions, ribbons are 95% credible intervals. Circles are training data, triangles are validation data.",
       x = "Year", y = "Total Coral Cover") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_output_dir, "validation_plot_no_pooling.pdf"), 
       plot = validation_plot, width = 15, height = 20, device = "pdf")

# Save the main validation plot
cat("Validation plot saved to:", file.path(plot_output_dir, "validation_plot_no_pooling.pdf"), "\n")