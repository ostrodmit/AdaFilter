# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The approach is to fit from the observations a time-invariant filter which reproduces the vector of observations with a small error, and at the same time has small norm in the DFT domain. This is done by solving a convex optimization problem.

The resulting estimator is adaptive to the unknown shift-invariant structure of the signal, encompassing, in particular, signals 
coming from low-dimensional shift-invariant linear subspaces (equivalently, discretized solutions of low-order ODEs and PDEs), and, more generally, signals well-approximated in such subspaces. In  particular, this includes the classical case of smooth estimation on the grid.

All the details concerning the theoretical motivation, statistical performance, and algorithmic implementation of our approach can be found in [these references][1,2,3]

[1] Ref1
[2] Ref2
[3] Ref3
