This folder contains code for computing an optimal control strategy for the COVID-19 outbreak in Ireland. The objective is to find a series of control measures which minimises the cost on the economic, while also ensuring that the hospital system is never overwhelmed. The code also for any number of switches in the level of disease control measures.

The file call_simulated_annealing.m runs the simulated annealing algorithm to find the global minimum of the objective function, as described in [1].
The function ode_solve_seir.m runs the SEIR model for a given control strategy, and computes the penalty function (economic cost with an additional large penalty for overwhelming the hospital system). 
The function ode_solve_seirX.m is used by the simulated annealing algorithm, which requires a function with only one output, the penalty function.

[1] O'Naraigh and Byrne, 'Piecewise-Constant Optimal Control Strategies for Controlling the Outbreak of COVID-19 in the Irish Population' (2020)
