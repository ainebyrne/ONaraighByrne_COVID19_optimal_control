This folder contains code for fitting a SEIR model to epidemiological data for the COVID-19 outbreak in Ireland.

The file call_simulated_annealing_parameters.m runs the simulated annealing algorithm to find the global minimum of the L-2 norm, as described in [1].
The function ode_solve_seir_parameters_betas.m runs the SEIR model for a given set of parameters, and computes the penalty function. 
The function ode_solve_seir_parameters_betasX.m is used by the simulated annealing algorithm, which requires a function with only one output, the penalty function.
The function calculate_R0.m computes the effective reproduction number for a given set of model parameters.
The epidemiological data, obtained from [2], is contained in CovidStatisticsProfileHPSCIrelandOpenData.csv

[1] O'Naraigh and Byrne, 'Piecewise-Constant Optimal Control Strategies for Controlling the Outbreak of COVID-19 in the Irish Population' (2020)
[2] Ordnance Survey Ireland. CovidStatisticsProfileHPSCIrelandOpenData. https://data.gov.ie/dataset/covidstatisticsprofilehpscirelandopendata. Visited: 15th June 2020.
