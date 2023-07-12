# Bioprocessing

[![Build Status](https://github.com/dfabianus/Bioprocessing.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/dfabianus/Bioprocessing.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Build Status](https://travis-ci.com/dfabianus/Bioprocessing.jl.svg?branch=master)](https://travis-ci.com/dfabianus/Bioprocessing.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/dfabianus/Bioprocessing.jl?svg=true)](https://ci.appveyor.com/project/dfabianus/Bioprocessing-jl)
[![Coverage](https://codecov.io/gh/dfabianus/Bioprocessing.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dfabianus/Bioprocessing.jl)
[![Coverage](https://coveralls.io/repos/github/dfabianus/Bioprocessing.jl/badge.svg?branch=master)](https://coveralls.io/github/dfabianus/Bioprocessing.jl?branch=master)


# Functionality of this package
- Bioprocess model library
    - mechanistic models
    - mass balances
    - elemental balances
    - reaction kinetics
    - reactor models and mass transfer
    - sensor models
    - data driven modeling
    - datatypes (AlgebraicEQs, ODE, PDE, DAE, SDE)
- Observer and controller library 
    - State estimation (KF, EKF, UKF, PF, MHE, ...)
    - PID, MPC, linear optimal control ...
    - data pre- and post-processing algorithms (Hampel, Smoothing, statistics)
    - Detection algorithms (Batch end, optimal point of harvest)
    - datatypes (ODE, discrete time DE)
- Model analysis and transformation 
    - sensitivity, identifiability, observability, jacobian
    - fisher information matrix
    - stiffness
    - controllability
    - Model reduction, inversion, linearization
    - datatypes (Linear algebra, nonlinear algorithms/functions, no measurements)
- Model simulation
    - simulating the ODE systems with ODE solvers
    - open loop plant (feedforward control, observer)
    - closed loop plant with (feedback control, observer)
    - bayesian stochastic simulations with parameter distributions
    - datatype (timeseries data)
- Model discovery and parameterization
    - Cost/Loss functions (NRMSE, ...)
    - optimization problem
    - bayesian maximum likelihood parameter estiamtion
    - SINDy, other symbolic regression algorithm
    - machine learning, regression, classification, DoE
    - model based design of experiments (mbDOE)

# Not in this package
- Software in the loop automation
- connection to external devices
- REST API connection to Lucullus