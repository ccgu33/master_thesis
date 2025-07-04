functions {
  vector coral_model(real t, vector y, vector theta_pars) {
    vector[2] dydt;
    real C = fmax(1e-6, fmin(1.0-1e-6, y[1]));
    real B = fmax(1e-6, fmin(1.0-1e-6, y[2]));

    real alpha = theta_pars[1];
    real beta = theta_pars[2];
    real gamma = theta_pars[3];
    real mu = theta_pars[4];

    real available_space = fmax(1e-6, 1.0 - (C + B));

    dydt[1] = alpha * C * available_space - beta * C + gamma * B;
    dydt[2] = beta * C - gamma * B - mu * B;

    return dydt;
  }
}

data {
  int<lower=1> N_years;
  int<lower=1> N_reefs;
  int<lower=1> N_sites;
  array[N_sites] int<lower=1, upper=N_reefs> reef_map;

  // Time setup
  real<lower=0> t0;
  array[N_years] real<lower=t0> ts;
  int<lower=1> no_t_gd;
  array[no_t_gd] real t_gd;
}

generated quantities {
  // Population level parameters
  vector[4] mu_pop;
  vector[4] tau_pop;
  cholesky_factor_corr[4] L_Omega_pop;
  
  // Reef level deviations
  matrix[4, N_reefs] z_reef;
  
  // Site level deviations
  vector[4] tau_site;
  cholesky_factor_corr[4] L_Omega_site;
  matrix[4, N_sites] z_site;

  // Observation error and initial conditions
  vector[N_sites] sigma;
  vector[N_sites] initial_cover;

  // Sample from uniform priors 
  for (p in 1:4) {
    mu_pop[p] = uniform_rng(0, 10);
    tau_pop[p] = uniform_rng(0, 5);
    tau_site[p] = uniform_rng(0, 5);
  }
  L_Omega_pop = lkj_corr_cholesky_rng(4, 1.0);
  L_Omega_site = lkj_corr_cholesky_rng(4, 1.0);

  // Sample random deviations
  for (r in 1:N_reefs) {
    for (p in 1:4) {
      z_reef[p, r] = normal_rng(0, 1);
    }
  }
  
  for (s in 1:N_sites) {
    for (p in 1:4) {
      z_site[p, s] = normal_rng(0, 1);
    }
    sigma[s] = uniform_rng(0, 3);
    initial_cover[s] = uniform_rng(0.01, 0.99);
  }

  // Build reef-level parameters
  matrix[4, N_reefs] mu_reef;
  matrix[4, N_reefs] temp_reef = rep_matrix(mu_pop, N_reefs) + diag_pre_multiply(tau_pop, L_Omega_pop) * z_reef;
  for (r in 1:N_reefs) {
    for (p in 1:4) {
      mu_reef[p, r] = fmax(1e-6, temp_reef[p, r]);
    }
  }

  // Build site-level parameters
  matrix[4, N_sites] theta;
  for (s in 1:N_sites) {
    vector[4] temp_site = mu_reef[, reef_map[s]] + diag_pre_multiply(tau_site, L_Omega_site) * z_site[,s];
    for (p in 1:4) {
      theta[p, s] = fmax(1e-6, temp_site[p]);
    }
  }

  // Set up initial conditions
  array[N_sites] vector[2] y0_sites;
  for (s in 1:N_sites) {
    real N0 = initial_cover[s];
    y0_sites[s, 1] = fmax(1e-6, fmin(1.0-1e-6, N0 * N0));
    y0_sites[s, 2] = fmax(1e-6, fmin(1.0-1e-6, N0 * (1.0 - N0)));
  }

  // Generate simulated data from priors
  array[N_years, N_sites] real Y_rep;
  for (s in 1:N_sites) {
    array[N_years] vector[2] y_hat = ode_bdf(coral_model, y0_sites[s], t0, ts, theta[, s]);
    for (t in 1:N_years) {
      real total_cover = fmax(1e-6, fmin(1.0-1e-6, y_hat[t, 1] + y_hat[t, 2]));
      Y_rep[t, s] = fmax(1e-6, fmin(1.0-1e-6, normal_rng(total_cover, sigma[s])));
    }
  }

  array[N_sites, no_t_gd] real y_smooth;
  if (no_t_gd > 0) {
    for (s in 1:N_sites) {
      array[no_t_gd] vector[2] y_traj = ode_bdf(coral_model, y0_sites[s], t0, t_gd, theta[, s]);
      for (t_grid in 1:no_t_gd) {
        y_smooth[s, t_grid] = fmax(1e-6, fmin(1.0-1e-6, y_traj[t_grid, 1] + y_traj[t_grid, 2]));
      }
    }
  }

  // Save parameters for checking
  vector[4] mu_pop_check = mu_pop;
  vector[4] tau_pop_check = tau_pop;
  vector[4] tau_site_check = tau_site;
  vector[N_sites] sigma_check = sigma;
  matrix[4, N_sites] theta_check = theta;
  vector[N_sites] initial_cover_check = initial_cover;
}
