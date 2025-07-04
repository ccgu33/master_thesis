functions {
  vector coral_model(real t, vector y, vector theta_pars) {
    vector[2] dydt;
    real C = fmax(0.0, fmin(1.0, y[1]));
    real B = fmax(0.0, fmin(1.0, y[2]));
    real alpha = theta_pars[1];
    real beta  = theta_pars[2];
    real gamma = theta_pars[3];
    real mu    = theta_pars[4];
    real available_space = fmax(0.0, 1.0 - (C + B));
    dydt[1] = alpha * C * available_space - beta * C + gamma * B;
    dydt[2] = beta * C - gamma * B - mu * B;
    return dydt;
  }
}

data {
  int<lower=1> N_years;
  int<lower=1> N_reefs;
  int<lower=1> N_sites;
  int<lower=0> N_obs_total;
  matrix[N_years, N_sites] Y_obs;
  matrix[N_years, N_sites] W_obs; // Forgetting factor weights
  int<lower=0, upper=1> Y_obs_present[N_years, N_sites];
  real<lower=0> t0;
  array[N_years] real<lower=t0> ts;
  array[N_sites] int<lower=1, upper=N_reefs> reef_map;
  real<lower=0> rtol;
  real<lower=0> atol;
  int<lower=1> max_num_steps;
  int<lower=1> no_t_gd;
  array[no_t_gd] real t_gd;
}

parameters {
  // Population level
  vector<lower=0>[4] mu_pop;     // Population means for alpha, beta, gamma, mu
  vector<lower=0>[4] tau_pop;    // How spread out reefs are (scale/width)
  cholesky_factor_corr[4] L_Omega_pop;  // Population correlation 
  
  // Reef level
  matrix[4, N_reefs] z_reef;     // Which specific direction each reef deviates (standardized units)
  
  // Site level
  vector<lower=0>[4] tau_site;   // Site-level standard deviations within reefs
  cholesky_factor_corr[4] L_Omega_site;  // Site-level correlation 
  matrix[4, N_sites] z_site;     // Raw site-level deviations
  
  // Observation level
  vector<lower=0>[N_sites] sigma;        // Observation noise 
  vector<lower=0.01, upper=0.99>[N_sites] y0_init_C;  // Initial coral cover proportion at each site
}

transformed parameters {
  // Reef level
  matrix<lower=0>[4, N_reefs] mu_reef;
  { matrix[4, N_reefs] mu_reef_raw = rep_matrix(mu_pop, N_reefs) + diag_pre_multiply(tau_pop, L_Omega_pop) * z_reef;
    for (r in 1:N_reefs) { for (p in 1:4) { mu_reef[p, r] = fmax(1e-6, mu_reef_raw[p, r]); } }
  }

  // Site level
  matrix<lower=0>[4, N_sites] theta;
   { matrix[4, N_sites] theta_raw;
    for (s in 1:N_sites) {
       theta_raw[,s] = mu_reef[, reef_map[s]] + diag_pre_multiply(tau_site, L_Omega_site) * z_site[,s];
       for (p in 1:4) { theta[p, s] = fmax(1e-6, theta_raw[p, s]); }
     }
   }

  // Initial conditions
  array[N_sites] vector[2] y0_sites;
  for (s in 1:N_sites) {
    real N0 = y0_init_C[s]; 
    y0_sites[s, 1] = fmax(1e-6, fmin(1.0 - 1e-6, N0 * N0));        // C(0)
    y0_sites[s, 2] = fmax(1e-6, fmin(1.0 - 1e-6, N0 * (1.0 - N0))); //B(0)
  }
}

model {
  // Population level priors: Student-t distributions selected via prior predictive checks
  mu_pop[1] ~ student_t(4, 1.6, 0.8) T[0, ];   // Alpha
  mu_pop[2] ~ student_t(4, 0.9, 0.5) T[0, ];   // Beta
  mu_pop[3] ~ student_t(4, 0.2, 0.2) T[0, ];   // Gamma 
  mu_pop[4] ~ student_t(4, 0.35, 0.25) T[0, ]; // Mu 
  tau_pop ~ student_t(4, 0, 0.5) T[0, ];       // Population-level variation 
  L_Omega_pop ~ lkj_corr_cholesky(2);          // Population-level correlation
  
  // Reef level priors
  to_vector(z_reef) ~ std_normal();
  
  // Site level priors
  tau_site ~ student_t(4, 0, 0.5) T[0, ];      // Site-level variation within reefs
  L_Omega_site ~ lkj_corr_cholesky(2);         // Site-level correlation between parameters
  to_vector(z_site) ~ std_normal();
  
  // Observation level priors
  sigma ~ student_t(4, 0, 0.2) T[0, ];         // Observation noise 
  y0_init_C ~ beta(2, 2);                      // Initial coral cover proportion (centered around 0.5)

  // Likelihood
  for (s in 1:N_sites) {
    array[N_years] vector[2] y_hat = ode_bdf_tol(coral_model, y0_sites[s], t0, ts, rtol, atol, max_num_steps, theta[, s]);
    for (t in 1:N_years) {
      if (Y_obs_present[t, s] == 1) {
        real pred_cov = fmax(0.0, fmin(1.0, y_hat[t, 1] + y_hat[t, 2])); // Total coral cover = C + B
        target += W_obs[t, s] * normal_lpdf(Y_obs[t, s] | pred_cov, sigma[s]);
      }
    }
  }
}

generated quantities {
  // Correlation matrices
  matrix[4, 4] Omega_pop = multiply_lower_tri_self_transpose(L_Omega_pop);   // Population-level parameter correlations
  matrix[4, 4] Omega_site = multiply_lower_tri_self_transpose(L_Omega_site); // Site-level parameter correlations
  
  // Posterior predictive checks, generate replicated datasets to assess model fit
  array[N_years, N_sites] real Y_rep;
  
  // Fitted trajectories, model predictions on dense time grid for visualization and analysis
  array[N_sites, no_t_gd] real y_fit_tot;

  for (s in 1:N_sites) {
     // Solve ODE model
     array[N_years] vector[2] y_hat_rep = ode_bdf_tol(coral_model, y0_sites[s], t0, ts, rtol, atol, max_num_steps, theta[, s]);
     
     // Generate replicated coral cover observations for posterior predictive checking
     for (t in 1:N_years) {
        if (Y_obs_present[t, s] == 1) {
          real pred_cov = fmax(0.0, fmin(1.0, y_hat_rep[t, 1] + y_hat_rep[t, 2])); 
          Y_rep[t, s] = fmax(0.0, fmin(1.0, normal_rng(pred_cov, sigma[s])));       
        } else {
          Y_rep[t, s] = -999;  // Missing data indicator
        }
     }
     
     if (no_t_gd > 0) {
       array[no_t_gd] vector[2] y_fit_comp = ode_bdf_tol(coral_model, y0_sites[s], t0, t_gd, rtol, atol, max_num_steps, theta[, s]);
       for (t_grid in 1:no_t_gd) {
         y_fit_tot[s, t_grid] = fmax(0.0, fmin(1.0, y_fit_comp[t_grid, 1] + y_fit_comp[t_grid, 2]));
       }
     }
  }
} 