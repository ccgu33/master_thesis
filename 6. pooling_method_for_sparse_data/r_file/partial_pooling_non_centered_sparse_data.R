library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)
library(stats) # For set.seed

inputdatafile_1 <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
stan_file_path <- "/Users/chengu/Documents/master thesis/final_final_final_code/stan_file_for_pooling_method/partial_pooling_non_centered_withoutforgetting_0703.stan"
rds_output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
rds_data_output_file <- file.path(rds_output_dir, "stan_data_partial_noncentered_sparse.rds")
rds_fit_output_file <- file.path(rds_output_dir, "stan_fit_partial_noncentered_sparse.rds")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(123) # For reproducible random sampling of reefs

# step 1: load and clean data
data_raw_full <- read.table(inputdatafile_1, header = TRUE, sep = ";", 
                           stringsAsFactors = FALSE, quote = "")

# select relevant columns and process data types
data_raw <- data_raw_full %>%
  select(REEF_NAME, SAMPLE_DATE, SITE_NO, GROUP_CODE, COVER, LATITUDE, LONGITUDE) %>%
  mutate(
    COVER = as.numeric(gsub("[^0-9.]", "", COVER)),
    SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%Y-%m-%d"),
    YEAR = as.numeric(format(SAMPLE_DATE, "%Y"))
  ) %>%
  filter(!is.na(YEAR), !is.na(COVER))

# calculate coral cover per site (proportion)
site_coral <- data_raw %>%
  filter(GROUP_CODE %in% c("Hard Coral")) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarize(
    site_coral = if(all(!is.na(COVER))) sum(COVER, na.rm = TRUE) / 100 else NA, 
    .groups = "drop"
  ) %>%
  mutate(site_coral = ifelse(is.finite(site_coral), site_coral, NA))

# step 2: select 15 reefs with complete data from 1995 to 2004
site_coral_filtered <- site_coral %>% 
  filter(YEAR >= 1995, YEAR <= 2004)

sites_per_reef_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(site_count = n_distinct(SITE_NO), .groups = "drop")

# find reefs with exactly 3 sites for all 10 years
reefs_with_complete_site_counts <- sites_per_reef_year %>% 
  group_by(REEF_NAME) %>%
  summarize(
    year_count = n_distinct(YEAR), 
    all_have_3_sites = all(site_count == 3), 
    complete_years = year_count == 10
  ) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

# get unique reef locations
complete_reef_locations <- data_raw_full %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% reefs_with_complete_site_counts) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

# step 3: select 15 reefs spread across latitude and longitude
# start with northernmost and southernmost reefs
extreme_reefs <- list(
  min_lat_reef = complete_reef_locations$REEF_NAME[which.min(complete_reef_locations$LATITUDE)],
  max_lat_reef = complete_reef_locations$REEF_NAME[which.max(complete_reef_locations$LATITUDE)]
)

selected_reefs_names <- unique(c(extreme_reefs$min_lat_reef, extreme_reefs$max_lat_reef))
target_n_reefs <- 15
reefs_needed <- target_n_reefs - length(selected_reefs_names)

# fill in remaining reefs by spreading across longitude
lat_extreme_longitudes <- complete_reef_locations$LONGITUDE[
  complete_reef_locations$REEF_NAME %in% selected_reefs_names
]
min_long <- min(lat_extreme_longitudes)
max_long <- max(lat_extreme_longitudes)
target_longitudes <- seq(min_long, max_long, length.out = reefs_needed + 2)[
  (1 + 1):(reefs_needed + 1)
]

remaining_reefs_df <- complete_reef_locations %>% 
  filter(!REEF_NAME %in% selected_reefs_names)

additional_reefs <- character(reefs_needed)
available_indices <- seq_len(nrow(remaining_reefs_df))

for (i in 1:reefs_needed) {
  distances <- abs(remaining_reefs_df$LONGITUDE[available_indices] - target_longitudes[i])
  closest_rel_idx <- which.min(distances)
  closest_abs_idx <- available_indices[closest_rel_idx]
  
  closest_reef_name <- remaining_reefs_df$REEF_NAME[closest_abs_idx]
  additional_reefs[i] <- closest_reef_name
  available_indices <- available_indices[-closest_rel_idx]
}

selected_reefs_names <- c(selected_reefs_names, additional_reefs[nzchar(additional_reefs)])
selected_reefs_names <- head(selected_reefs_names, target_n_reefs)



# step 4: map reef names to numbers for stan
reef_id_to_name_map_initial <- complete_reef_locations %>%
  filter(REEF_NAME %in% selected_reefs_names) %>%
  arrange(LONGITUDE) %>%
  mutate(Reef_ID = 1:n()) %>%
  select(Reef_ID, REEF_NAME)

# step 5: prepare initial data with all sites
all_sites_data <- site_coral %>%
  filter(REEF_NAME %in% selected_reefs_names, YEAR >= 1995, YEAR <= 2004) %>%
  left_join(reef_id_to_name_map_initial, by = "REEF_NAME") %>%
  mutate(SITE_NO = as.character(SITE_NO)) %>%
  arrange(Reef_ID, SITE_NO, YEAR) %>%
  mutate(initial_site_id = match(paste(Reef_ID, SITE_NO, sep="_"), 
                                unique(paste(Reef_ID, SITE_NO, sep="_"))))

# step 6: create sparse data structure (5x3, 5x2, 5x1 sites per reef)
# randomly assign reefs to have 3, 2, or 1 site(s)
unique_selected_reefs <- reef_id_to_name_map_initial$REEF_NAME
reefs_3_sites <- sample(unique_selected_reefs, 5)
remaining_reefs_1 <- setdiff(unique_selected_reefs, reefs_3_sites)
reefs_2_sites <- sample(remaining_reefs_1, 5)
reefs_1_site <- setdiff(remaining_reefs_1, reefs_2_sites)



# filter the data to create sparse structure
reduced_sites_data <- all_sites_data %>%
  filter(
    # Keep all sites for reefs_3_sites 
    # Keep only sites '1' and '2' for reefs_2_sites  
    !(REEF_NAME %in% reefs_2_sites & SITE_NO == "3") &
    !(REEF_NAME %in% reefs_1_site & (SITE_NO == "2" | SITE_NO == "3"))
  )

# step 7: prepare final data with updated site ids
site_id_mapping <- reduced_sites_data %>%
  distinct(Reef_ID, SITE_NO) %>%
  arrange(Reef_ID, SITE_NO) %>%
  mutate(final_site_id = 1:n())

final_data <- reduced_sites_data %>%
  left_join(site_id_mapping, by = c("Reef_ID", "SITE_NO"))

# recalculate dimensions based on final data
N_sites <- max(final_data$final_site_id)
N_reefs <- max(final_data$Reef_ID)
unique_years <- sort(unique(final_data$YEAR))
years_relative <- unique_years - min(unique_years)
N_years <- length(years_relative)

t0 <- 0
ts <- years_relative + 1e-9  # tiny number to avoid ode solver issues
# create final reef map for stan
reef_map_stan <- site_id_mapping %>%
  arrange(final_site_id) %>%
  pull(Reef_ID)

# step 8: transform to matrix format for stan (with missing data handling)
coral_matrix <- final_data %>%
  select(YEAR, final_site_id, site_coral) %>%
  mutate(Year_Index = match(YEAR, unique_years)) %>%
  pivot_wider(
    id_cols = Year_Index, 
    names_from = final_site_id, 
    values_from = site_coral, 
    names_prefix = "site_"
  ) %>%
  arrange(Year_Index) %>% 
  select(-Year_Index) %>% 
  as.matrix()

# Create indicator matrix for missing data
obs_indicator <- ifelse(is.na(coral_matrix), 0, 1)
N_obs_total <- sum(obs_indicator)



# Replace NAs with placeholder for Stan
coral_matrix[is.na(coral_matrix)] <- 0

t_gd <- seq(from = min(years_relative), to = max(years_relative), length.out = 100)
no_t_gd <- length(t_gd)
t_gd <- pmax(t_gd + 1e-9, ts[1])

# Step 9: Assemble Stan data list
stan_data_list <- list(
  N_years = N_years,
  N_reefs = N_reefs,
  N_sites = N_sites,
  N_obs_total = N_obs_total,
  Y_obs = coral_matrix,
  Y_obs_present = obs_indicator,
  t0 = t0,
  ts = ts,
  reef_map = reef_map_stan,
  rtol = 1e-6,
  atol = 1e-6,
  max_num_steps = 100000,
  no_t_gd = no_t_gd,
  t_gd = t_gd
)

saveRDS(stan_data_list, rds_data_output_file)

# Save reef mapping for plotting
rds_mapping_output_file <- file.path(rds_output_dir, "reef_id_mapping_sparse.rds")
saveRDS(
  reef_id_to_name_map_initial %>% 
  left_join(complete_reef_locations, by="REEF_NAME") %>%
  select(Reef_ID, REEF_NAME, LATITUDE, LONGITUDE),
  rds_mapping_output_file
)

compiled_model <- stan_model(stan_file_path)

# Step 11: Run MCMC sampling and save results
fit <- sampling(
  compiled_model,
  data = stan_data_list,
  seed = 123,
  chains = 4,
  iter = 1000,
  warmup = 500
 
)

saveRDS(fit, rds_fit_output_file)

