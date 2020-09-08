This repository contains code for computing an optimal control strategy for the COVID-19 outbreak in Ireland. The objective is to find a series of control measures which minimises the cost on the economic, while also ensuring that the hospital system is never overwhelmed. The corresponding results are presented in [1]. 

The folder Fitting-SEIR-model-to-data contains code for fitting the SEIR model to the epidemiological data, as described in Section 2 of [1].
The folder Optimal-control-with-one-switch contains code for determining the optimal control strategy with one switch in the level of control measures. The optimal strategy for this scenarios is described in Section 4.1 of [1].
The code in Optimal-control-with-n-switches contains allows for multiple switches in the level of control measures. The results for this control strategy are presented in Section 4.2 of [1].

[1] O'Naraigh and Byrne, 'Piecewise-Constant Optimal Control Strategies for Controlling the Outbreak of COVID-19 in the Irish Population' (2020)
