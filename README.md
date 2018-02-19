# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The approach is to fit a (data-dependent and hence non-linear) time-invariant filter which reproduces the vector of observations with a small error, and at the same time has a small l1-norm in the DFT domain. This is done by solving an appropriate convex optimization problem.

The resulting estimator is adaptive to the unknown shift-invariant structure of the signal, encompassing, in particular, signals 
coming from low-dimensional shift-invariant linear subspaces (equivalently, discretized solutions of low-order ODEs and PDEs), and, more generally, signals well-approximated in such subspaces. In  particular, this includes the classical case of smooth estimation on the grid.

All the details concerning the theore
