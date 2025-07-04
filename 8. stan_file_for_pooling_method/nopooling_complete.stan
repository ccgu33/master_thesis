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
  int<lower=1> N_sites;
  matrix[N_years, N_sites] Y_obs;
  int<lower=0, upper=1> Y_obs_present[N_years, N_sites];
  real<lower=0> t0;
  array[N_years] real<lower=t0> ts;
  real<lower=0> rtol;
  real<lower=0> atol;
  int<lower=1> max_num_steps;
  int<lower=1> no_t_gd;
  array[no_t_gd] real t_gd;
}

parameters {
  // Site level
  matrix<lower=0>[4, N_sites] theta;  
  
  // Observation level
  vector<lower=0>[N_sites] sigma;        // Observation noise 
  vector<lower=0.01, upper=0.99>[N_sites] y0_init_C;  // Initial coral cover proportion at each site
}

transformed parameters {
  // Initial conditions
  array[N_sites] vector[2] y0_sites;
  for (s in 1:N_sites) {
    real N0 = y0_init_C[s]; 
    y0_sites[s, 1] = fmax(1e-6, fmin(1.0 - 1e-6, N0 * N0));        // C(0)
    y0_sites[s, 2] = fmax(1e-6, fmin(1.0 - 1e-6, N0 * (1.0 - N0))); //B(0)
  }
}

model {
  // Site level priors
  for (s in 1:N_sites) {
    theta[1, s] ~ student_t(4, 1.6, 0.8) T[0, ];   // Alpha
    theta[2, s] ~ student_t(4, 0.9, 0.5) T[0, ];   // Beta
    theta[3, s] ~ student_t(4, 0.2, 0.2) T[0, ];   // Gamma 
    theta[4, s] ~ student_t(4, 0.35, 0.25) T[0, ]; // Mu 
  }

  // Observation level priors
  sigma ~ student_t(4, 0, 0.1) T[0, ];         // Observation noise 
  y0_init_C ~ beta(2, 2);                      // Initial coral cover proportion (centered around 0.5)

  // Likelihood
  for (s in 1:N_sites) {
    array[N_years] vector[2] y_hat = ode_bdf_tol(coral_model, y0_sites[s], t0, ts, rtol, atol, max_num_steps, theta[, s]);
    for (t in 1:N_years) {
      if (Y_obs_present[t, s] == 1) {
         real predicted_total_cover = y_hat[t, 1] + y_hat[t, 2]; // Total coral cover = C + B
         predicted_total_cover = fmax(0.0, fmin(1.0, predicted_total_cover));
         Y_obs[t, s] ~ normal(predicted_total_cover, sigma[s]);
      }
    }
  }
}

generated quantities {
  array[N_years, N_sites] real Y_rep;
  
  array[N_sites, no_t_gd] real y_fit_tot;

  for (s in 1:N_sites) {
     // Solve ODE model
     array[N_years] vector[2] y_hat_rep = ode_bdf_tol(coral_model, y0_sites[s], t0, ts, rtol, atol, max_num_steps, theta[, s]);
     for (t in 1:N_years) {
        if (Y_obs_present[t, s] == 1) {
          real pred_cov = fmax(0.0, fmin(1.0, y_hat_rep[t, 1] + y_hat_rep[t, 2])); 
          Y_rep[t, s] = fmax(0.0, fmin(1.0, normal_rng(pred_cov, sigma[s])));      
        } else {
          Y_rep[t, s] = -999;  
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
