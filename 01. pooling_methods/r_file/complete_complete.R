library(rstan)    
library(dplyr)   
library(tidyr)    
library(magrittr)

inputdatafile_1 <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
stan_file_path <- "/Users/chengu/Documents/master thesis/final_final_final_code/stan_file_for_pooling_method/complete_pooling_complete_0703.stan"
rds_output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
rds_data_output_file <- file.path(rds_output_dir, "stan_data_complete_pooling_noraml_2.rds")
rds_fit_output_file <- file.path(rds_output_dir, "stan_fit_complete_pooling_normal_2.rds")

options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)             

# Step 1: Load and clean data
data_raw_full <- read.table(
  inputdatafile_1, 
  header = TRUE, 
  sep = ";", 
  stringsAsFactors = FALSE, 
  quote = ""
) 

# Select relevant columns and process data types
data_raw <- data_raw_full %>%
  select(
    REEF_NAME, SAMPLE_DATE, SITE_NO, 
    GROUP_CODE, COVER, LATITUDE, LONGITUDE
  ) %>% 
  mutate(
    COVER = as.numeric(gsub("[^0-9.]", "", COVER)),
    SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%Y-%m-%d"),
    YEAR = as.numeric(format(SAMPLE_DATE, "%Y"))
  )

# Calculate coral cover per site
site_coral <- data_raw %>%
  # Filter for Hard Coral data only
  filter(GROUP_CODE %in% c("Hard Coral")) %>% 
  # Group by reef, year, and site
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  # Calculate total coral cover (divide by 100 to get proportion)
  summarize(site_coral = sum(COVER) / 100, .groups = "drop")

# Step 2: Select 15 reefs with complete data 1995-2004
site_coral_filtered <- site_coral %>% 
  filter(YEAR >= 1995, YEAR <= 2004)

# Count number of sites per reef per year
sites_per_reef_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(
    site_count = n_distinct(SITE_NO), 
    .groups = "drop"
  )

# reefs with exactly 3 sites for all 10 years
reefs_with_complete_site_counts <- sites_per_reef_year %>% 
  group_by(REEF_NAME) %>%
  summarize(
    year_count = n_distinct(YEAR),
    all_have_3_sites = all(site_count == 3),
    complete_years = year_count == 10
  ) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

complete_reefs_final <- reefs_with_complete_site_counts

if(length(complete_reefs_final) < 15) {
  stop("Fewer than 15 complete reefs available.")
}

complete_reef_locations <- data_raw_full %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% complete_reefs_final) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

# Start with northernmost and southernmost reefs
extreme_reefs <- list(
  min_lat_reef = complete_reef_locations$REEF_NAME[which.min(complete_reef_locations$LATITUDE)],
  max_lat_reef = complete_reef_locations$REEF_NAME[which.max(complete_reef_locations$LATITUDE)]
)

selected_reefs_names <- unique(c(extreme_reefs$min_lat_reef, extreme_reefs$max_lat_reef))
target_n_reefs <- 15
reefs_needed <- target_n_reefs - length(selected_reefs_names)

# Fill in remaining reefs by spreading across longitude
if (reefs_needed > 0) {
  lat_extreme_longitudes <- complete_reef_locations$LONGITUDE[
    complete_reef_locations$REEF_NAME %in% selected_reefs_names
  ]
  min_long <- min(lat_extreme_longitudes)
  max_long <- max(lat_extreme_longitudes)
  #skip the first and last 
  target_longitudes <- seq(min_long, max_long, length.out = reefs_needed + 2)[
    (1 + 1):(reefs_needed + 1)
  ]
  
  remaining_reefs_df <- complete_reef_locations %>% 
    filter(!REEF_NAME %in% selected_reefs_names)
  
  additional_reefs <- character(reefs_needed)
  
  for (i in 1:reefs_needed) {
    if(nrow(remaining_reefs_df) == 0) break
    
    distances <- abs(remaining_reefs_df$LONGITUDE - target_longitudes[i])
    closest_idx <- which.min(distances)
    closest_reef_name <- remaining_reefs_df$REEF_NAME[closest_idx]
    
    additional_reefs[i] <- closest_reef_name
    remaining_reefs_df <- remaining_reefs_df[-closest_idx, ]
  }
  
  selected_reefs_names <- c(selected_reefs_names, additional_reefs[nzchar(additional_reefs)])
}

selected_reefs_names <- head(unique(selected_reefs_names), target_n_reefs)

if(length(selected_reefs_names) < target_n_reefs) {
  stop("Could not select 15 unique reefs.")
}

# Map reef names to numbers for Stan
reef_id_to_name_map <- complete_reef_locations %>%
  filter(REEF_NAME %in% selected_reefs_names) %>%
  arrange(LONGITUDE) %>%
  mutate(Reef_ID = 1:n()) %>%
  select(Reef_ID, REEF_NAME)

# Step 3: Prepare data for Stan
hierarchical_data_long <- site_coral %>%
  filter(REEF_NAME %in% selected_reefs_names, YEAR >= 1995, YEAR <= 2004) %>%
  left_join(reef_id_to_name_map, by = "REEF_NAME") %>%
  arrange(Reef_ID, SITE_NO, YEAR) %>%
  group_by(REEF_NAME, SITE_NO) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup()

N_sites <- max(hierarchical_data_long$site_id)
unique_years <- sort(unique(hierarchical_data_long$YEAR))
years_relative <- unique_years - min(unique_years)
N_years <- length(years_relative)

t0 <- 0
ts_raw <- years_relative
ts <- ts_raw + 1e-9  # Tiny number to avoid ODE solver issues

cat("N_years =", N_years, ", N_sites =", N_sites, "\n")

# Transform to matrix format for stan
Y_data_matrix <- hierarchical_data_long %>%
  select(YEAR, site_id, site_coral) %>%
  mutate(Year_Index = match(YEAR, unique_years)) %>%
  pivot_wider(
    id_cols = Year_Index, 
    names_from = site_id, 
    values_from = site_coral, 
    names_prefix = "site_"
  ) %>%
  arrange(Year_Index) %>% 
  select(-Year_Index) %>% 
  as.matrix()

Y_obs_present <- matrix(1, nrow = nrow(Y_data_matrix), ncol = ncol(Y_data_matrix))
N_obs_total <- sum(Y_obs_present)

# Time points for predictions
t_gd <- seq(from = min(ts_raw), to = max(ts_raw), length.out = 100)
no_t_gd <- length(t_gd)

if (no_t_gd > 0 && t_gd[1] <= t0) {
  t_gd <- t_gd + 1e-9
  t_gd <- pmax(t_gd, ts[1])
}

stan_data_list <- list(
  N_years = N_years,
  N_sites = N_sites,
  Y_obs = Y_data_matrix,
  Y_obs_present = Y_obs_present,
  t0 = t0,
  ts = ts,
  no_t_gd = no_t_gd,
  t_gd = t_gd
)



saveRDS(stan_data_list, rds_data_output_file)

cat("Compiling Stan model:", stan_file_path, "\n")
if (!file.exists(stan_file_path)) {
  stop("Stan file not found: ", stan_file_path)
}
compiled_model <- stan_model(stan_file_path)

# Step 5: Run MCMC sampling and save results
fit <- sampling(
  compiled_model,
  data = stan_data_list,
  chains = 4,
  iter = 2000,
  warmup = 1000
)

saveRDS(fit, "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds/stan_fit_complete_pooling.rds")