# SPADE_prediction

This repository contains R code and instructions to reproduce the results presented in the manuscript "A Mixed-Effect Dirichlet Multinomial Model for Evaluating Hospital Diagnostic Performance." The analysis consists of two main parts: Simulation and Real Data Analysis.

## Simulation Study
The script `Simulation.R` includes generation of the simulated dataset, maximum likelihood estimation, optimization methods, evaluation of results, including confidence and prediction intervals. `sim_results` generates the results shown in Table 1 and Figure 1 in the manuscript.

## Real Data Analysis
`Syntax_data_generation.R`: Method of generating a synthetic dataset.
`syn_cov.csv`: Contains the simulated covariate data.
`syn_ni_list`: Simulates the number of patient returns for three patient types.
`MLEstimation.R`: Performs maximum likelihood estimation.
`da_result.R`: Analyzes results and generates Table 2 and Figure 2 in the manuscript.
