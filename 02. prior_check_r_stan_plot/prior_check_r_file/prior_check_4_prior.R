library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(reshape2)

base_dir <- "/Users/chengu/Documents/master thesis/final_final_stan_and r_code"
inputdatafile_1 <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
stan_model_dir <- file.path(base_dir, "prior_check_stan_file")
rds_base_dir <- file.path(base_dir, "rds")
output_base_dir <- file.path(base_dir, "prior_check_output")

output_dir <- file.path(output_base_dir, "prior_run")
dir.create(output_dir, recursive = TRUE)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
n_prior_samples <- 1000
ode_time_epsilon <- 1e-9

#step1: load data
data_raw_full <- read.table(inputdatafile_1, header = TRUE, sep = ";", stringsAsFactors = FALSE, quote = "")
data_raw <- data_raw_full %>%
  select(REEF_NAME, SAMPLE_DATE, SITE_NO, GROUP_CODE, COVER, LATITUDE, LONGITUDE) %>%
  mutate(COVER = as.numeric(gsub("[^0-9.]", "", COVER)),
         SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%Y-%m-%d"),
         YEAR = as.numeric(format(SAMPLE_DATE, "%Y")))

site_coral <- data_raw %>%
  filter(GROUP_CODE %in% c("Hard Coral")) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarize(site_coral = sum(COVER, na.rm = FALSE) / 100, .groups = "drop") %>%
  mutate(site_coral = ifelse(is.finite(site_coral), site_coral, NA))

site_coral_filtered <- site_coral %>% filter(YEAR >= 1995, YEAR <= 2004)

sites_per_reef_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(site_count = n_distinct(SITE_NO), .groups = "drop")

reefs_with_complete_site_counts <- sites_per_reef_year %>% 
  group_by(REEF_NAME) %>%
  summarize(year_count = n_distinct(YEAR), 
            all_have_3_sites = all(site_count == 3), 
            complete_years = year_count == 10) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

reefs_with_na_values <- site_coral_filtered %>% 
  filter(REEF_NAME %in% reefs_with_complete_site_counts) %>%
  group_by(REEF_NAME) %>% 
  summarize(has_na = any(is.na(site_coral))) %>% 
  filter(has_na) %>% 
  pull(REEF_NAME)

complete_reefs_final <- setdiff(reefs_with_complete_site_counts, reefs_with_na_values)

complete_reef_locations <- data_raw_full %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% complete_reefs_final) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

extreme_reefs <- list(
  min_lat_reef = complete_reef_locations$REEF_NAME[which.min(complete_reef_locations$LATITUDE)],
  max_lat_reef = complete_reef_locations$REEF_NAME[which.max(complete_reef_locations$LATITUDE)]
)
selected_reefs_names <- unique(c(extreme_reefs$min_lat_reef, extreme_reefs$max_lat_reef))
target_n_reefs <- 15
reefs_needed <- target_n_reefs - length(selected_reefs_names)

longitudes <- complete_reef_locations$LONGITUDE[complete_reef_locations$REEF_NAME %in% selected_reefs_names]
target_longitudes <- seq(min(longitudes), max(longitudes), length.out = reefs_needed + 2)[2:(reefs_needed + 1)]
remaining_reefs_df <- complete_reef_locations %>% filter(!REEF_NAME %in% selected_reefs_names)
additional_reefs <- character(reefs_needed)

for (i in 1:reefs_needed) {
  distances <- abs(remaining_reefs_df$LONGITUDE - target_longitudes[i])
  closest_idx <- which.min(distances)
  closest_reef_name <- remaining_reefs_df$REEF_NAME[closest_idx]
  additional_reefs[i] <- closest_reef_name
  remaining_reefs_df <- remaining_reefs_df[-closest_idx, ]
}
selected_reefs_names <- c(selected_reefs_names, additional_reefs[nzchar(additional_reefs)])
selected_reefs_names <- head(unique(selected_reefs_names), target_n_reefs)

reef_id_to_name_map <- complete_reef_locations %>%
  filter(REEF_NAME %in% selected_reefs_names) %>%
  arrange(LONGITUDE) %>%
  mutate(Reef_ID = 1:n()) %>%
  select(Reef_ID, REEF_NAME)

print(reef_id_to_name_map)

# step 2: create stan data structure
N_reefs <- 15
N_sites <- 45  # 15 reefs × 3 sites per reef
N_years <- 10  # 1995-2004
reef_map <- rep(1:15, each = 3)  # Sites 1-3 = Reef 1, Sites 4-6 = Reef 2, etc.

# Time points
t0 <- 1995
ts <- c(1995:2004) + ode_time_epsilon
t_gd <- seq(1995, 2004, length.out = 50) + ode_time_epsilon

stan_data_structure <- list(
  N_years,
  N_reefs,
  N_sites,
  reef_map,
  t0,
  ts,
  max_num_steps = 1000,
  no_t_gd = length(t_gd),
  t_gd
)

#step 3: run prior checks
prior_model_files <- list(
  `Uniform` = file.path(stan_model_dir, "prior_check_flat.stan"),
  `Vague Normal` = file.path(stan_model_dir, "prior_check_vague_large_normal.stan"),
  `Tight Normal` = file.path(stan_model_dir, "prior_check_small_normal.stan"),
  `Student-t` = file.path(stan_model_dir, "prior_check_generic_weakly_informatice_student_t.stan")
)

all_prior_samples <- list()

for (model_name in names(prior_model_files)) {
  stan_file <- prior_model_files[[model_name]]
  compiled_model <- stan_model(stan_file)
  
  fit_prior <- sampling(compiled_model, 
                       data = stan_data_structure, 
                       chains = 1, 
                       iter = n_prior_samples, 
                       warmup = 0,
                       seed = 123 + match(model_name, names(prior_model_files)), 
                       algorithm = "Fixed_param")
  
  saveRDS(fit_prior, file = file.path(output_dir, paste0("prior_fit_", gsub(" ", "_", model_name), ".rds")))
  
  prior_samples <- rstan::extract(fit_prior)
  all_prior_samples[[model_name]] <- prior_samples
}

#step 4: create data frame for plotting
#a. create site_to_reef_mapping
site_to_reef_mapping <- data.frame(
  Site_ID = 1:stan_data_structure$N_sites,  # 1-45
  Reef_ID = stan_data_structure$reef_map  # 1-15
) %>%
  merge(reef_id_to_name_map, by = "Reef_ID")

#b. create combined_y_rep_all list to store data from all models
combined_y_rep_all <- list()
for (model_name in names(all_prior_samples)) { 
  y_rep_array <- all_prior_samples[[model_name]]$Y_rep # extract 3d array 1000 x 10 x 45

 #c. convert 3d array to long format
  y_rep_long <- reshape2::melt(y_rep_array, # 1000 x 10 x 45
                              varnames = c("Sim", "Year_Index", "Site_ID"), # Sim: 1-1000, Year_Index: 1-10, Site_ID: 1-45
                              value.name = "Cover")
  
  y_rep_long$Site_ID <- as.integer(y_rep_long$Site_ID)
  
  #d. merge with site_to_reef_mapping to get reef names
  y_rep_long_with_reef <- merge(y_rep_long, site_to_reef_mapping, by = "Site_ID")
  #e. add model name column
  y_rep_long_with_reef$Model <- model_name
  #f. convert year index to actual year
  y_rep_long_with_reef$Year <- ts[y_rep_long_with_reef$Year_Index]
  #g. store processed data for this model
  combined_y_rep_all[[model_name]] <- y_rep_long_with_reef
}
#h. combine all model data frames into one
combined_y_rep_df <- do.call(rbind, combined_y_rep_all)

#step 5: create comparison plots

# Plot 1: Overall density
overall_density_plot <- ggplot(combined_y_rep_df, aes(x = Cover, color = Model)) +
  geom_density(alpha = 0.7, linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Prior Predictive Distributions - Overall",
       subtitle = paste(n_prior_samples, "samples ×", stan_data_structure$N_years, "years ×", stan_data_structure$N_sites, "sites"),
       x = "Simulated Coral Cover", 
       y = "Density") +
  scale_color_brewer(palette = "Set1") +  
  theme_bw()

ggsave(file.path(output_dir, "compare_prior_pred_final_cover_density.png"), 
       plot = overall_density_plot, width = 9, height = 6)

# Plot 2: Overall histogram
overall_hist_plot <- ggplot(combined_y_rep_df, aes(x = Cover)) +
  geom_histogram(bins = 50, fill = "gray", alpha = 0.7) +
  facet_wrap(~ Model) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Prior Predictive Distributions - Overall",
       x = "Simulated Coral Cover", 
       y = "Count") +
  theme_bw()

ggsave(file.path(output_dir, "compare_prior_pred_final_cover_hist.png"), 
       plot = overall_hist_plot, width = 10, height = 7)

# Plot 3: Reef-specific density
reef_density_plot <- ggplot(combined_y_rep_df, aes(x = Cover, color = Model)) +
  geom_density(alpha = 0.8, linewidth = 0.7) +
  facet_wrap(~ REEF_NAME, ncol = 5) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Prior Predictive Distributions by Reef (All 15 Reefs)",
       subtitle = paste(n_prior_samples, "samples ×", stan_data_structure$N_years, "years × 3 sites per reef"),
       x = "Simulated Coral Cover", 
       y = "Density") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 8) +  
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 6)) 

ggsave(file.path(output_dir, "compare_prior_pred_reef_cover_density.png"), 
       plot = reef_density_plot, width = 18, height = 12) 