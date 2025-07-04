# Bayesian Inference of Coral Bleaching Dynamics

This repository contains the code and analysis for my master thesis on Bayesian inference of coral bleaching dynamics. The project uses hierarchical Bayesian models to analyze coral cover data and understand bleaching patterns across different reefs.

## Project Structure

The repository is organized as follows:

- **1. pooling_methods/**: Implementation of different pooling approaches
  - **r_file/**: R scripts for complete, no, and partial pooling models
  - **stan_file/**: Stan model files for Bayesian inference

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

## Usage

Each R script in the `r_file` directory corresponds to a different pooling method and can be run independently.

## Author

Chen Gu 