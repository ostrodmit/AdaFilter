# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The goal is to estimate a discrete-time signal 

The approach fits a data-dependent time-invariant filter by solving an appropriate convex optimization problem (well-structured SOCP).
The resulting estimator is adaptive to the unknown shift-invariant structure of the signal, encompassing, in particular, signals satisfying low-order ODEs and PDEs, and more generally, signals 

All the details concerning the theore
