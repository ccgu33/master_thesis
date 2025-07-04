# Bayesian Inference of Coral Bleaching Dynamics

This repository contains the code and analysis for my master thesis on Bayesian inference of coral bleaching dynamics. The project uses hierarchical Bayesian models to analyze coral cover data and understand bleaching patterns across different reefs.

## Repository Statistics

The repository contains a comprehensive collection of files for Bayesian analysis:
- **Total Files**: 95 files
- **Directories**: 11 main directories with multiple subdirectories
- **R Scripts**: 14 R files for data processing, analysis, and visualization
- **Stan Models**: 16 Stan model files for Bayesian inference
- **RDS Files**: 27 RDS files containing saved data and model fits

## Project Structure

The repository is organized into 11 main directories:

1. **1. pooling_methods/**: Implementation of different pooling approaches
   - **r_file/**: 5 R scripts for complete, no, and partial pooling models
   - **stan_file/**: 6 Stan model files for Bayesian inference
   - **rds/**: 14 RDS files with model data and fitted objects

2. **2. prior_check_r_stan_plot/**: Prior distribution analysis
   - **prior_check_r_file/**: R scripts for prior checks
   - **prior_check_stan_file/**: 4 Stan models with different prior specifications
   - **prior_check_rds/**: RDS files for prior check results

3. **3. posterier_predictive_check/**: Model validation
   - 3 R scripts for Bayesian p-values and posterior predictive checks

4. **4. leave_one_out_validation/**: Cross-validation analysis
   - **r file/**: 2 R scripts for LOO validation and plotting
   - **rds/**: RDS files with LOO comparison results

5. **5. mcmc_diagnostics/**: MCMC convergence analysis
   - R scripts for diagnosing MCMC sampling quality

6. **6. pooling_method_for_sparse_data/**: Sparse data analysis
   - **r_file/**: R script for partial pooling with sparse data
   - **rds/**: RDS files with sparse data model fits

7. **7. prediction_for_unseen_data/**: Out-of-sample prediction
   - **r_file/**: R script for predicting on new data
   - **rds/**: RDS files with prediction results

8. **8. stan_file_for_pooling_method/**: 6 Stan model files (duplicates from directory 1)

9. **9. plots/**: Visualization outputs
   - Subdirectories for different types of plots (prior checks, posterior predictions, etc.)

10. **10. master_thesis/**: Thesis document files

11. **11. data/**: Raw and processed data files

## Models

The project implements several Bayesian modeling approaches:
- Complete pooling: Assumes all reefs share identical parameters
- No pooling: Treats each reef as completely independent
- Partial pooling: Hierarchical model that balances between complete and no pooling
  - Centered parameterization
  - Non-centered parameterization
  - Sparse data handling

## Requirements

- R (>= 4.0.0)
- RStan
- dplyr
- tidyr
- magrittr
- loo (for leave-one-out cross-validation)
- bayesplot (for MCMC diagnostics and visualization)

## Usage

Each R script in the `r_file` directory corresponds to a different pooling method and can be run independently. The workflow typically involves:

1. Data preprocessing
2. Model fitting with Stan
3. MCMC diagnostics
4. Posterior predictive checks
5. Model comparison via LOO-CV
6. Visualization of results

## Author

Chen Gu
