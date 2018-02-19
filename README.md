# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The approach is to fit from the observations a time-invariant filter which *reproduces the vector of observations with a small error, and at the same time has small norm in the DFT domain.* This is done by solving a convex optimization problem (namely a well-structured SOCP).

The resulting estimator is adaptive to the unknown shift-invariant structure of the signal. Thus, we address the case of signals 
which belong to low-dimensional shift-invariant linear subspaces (equivalently, discretized solutions of low-order ODEs and PDEs), and, more generally, signals that are well-approximated in such subspaces. This includes the classical case of smooth estimation on the grid.

All the details regarding the theoretical motivation, statistical performance, and algorithmic implementation of the outlined approach can be found in [1,2,3]. We also plan to release the code reproducing the experimental results of those papers.

## Usage

## Documentation

## References

[1] Ref1
[2] Ref2
[3] Ref3
