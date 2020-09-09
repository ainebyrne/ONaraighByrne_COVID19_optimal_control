This folder contains code for fitting a SEIR model to epidemiological data for the COVID-19 outbreak in Ireland. Bootstrapping is used to create synthetic data and compute the confidence intervals for the best fitting parameters

The file generate_confidence_interval.m runs a simulated annealing algorithm (as described in [1]) to find the best fitting parameters for the actual data [2]. It then generated synthetic data to fit the model to and computes the confidence intervals for the best fitting parameter values.
The function ode_solve_seir_parameters_betas.m runs the SEIR model for a given set of parameters, and computes the penalty function. 
The function ode_solve_seir_parameters_betasX.m is used by the simulated annealing algorithm, which requires a function with only one output, the penalty function.
The function ode_solve_synthetic.m and ode_solve_synthetic_wrapper.m are the equivalent functions for the synthetic data. 
The function calculate_R0.m computes the effective reproduction number for a given set of model parameters.
The epidemiological data, obtained from [2], is contained in CovidStatisticsProfileHPSCIrelandOpenData.csv
The computed confidence intervals are reported in CIs.xlsx

[1] O'Naraigh and Byrne, 'Piecewise-Constant Optimal Control Strategies for Controlling the Outbreak of COVID-19 in the Irish Population' (2020)
[2] Ordnance Survey Ireland. CovidStatisticsProfileHPSCIrelandOpenData. https://data.gov.ie/dataset/covidstatisticsprofilehpscirelandopendata. Visited: 15th June 2020.
