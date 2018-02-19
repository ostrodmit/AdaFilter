# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The goal is to estimate a discrete-time signal on a finite regular grid.

The approach fits a data-dependent time-invariant filter by solving an appropriate convex optimization problem.
The resulting estimator is adaptive to the unknown shift-invariant structure of the signal, encompassing, in particular, signals 
coming from low-dimensional shift-invariant linear subspaces (equivalently, discretized solutions of low-order ODEs and PDEs), and, more generally, signals well-approximated in such subspaces. In  particular, this includes the classical case of smooth estimation on the grid.

All the details concerning the theore
