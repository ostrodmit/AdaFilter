# AdaFilter — Signal Denoising via Adaptive Filtering

AdaFilter provides a **Matlab** implementation of struture-adaptive denoising of discrete-time signals and images.

## Short description

The approach is to fit from the observations a time-invariant filter which *reproduces the vector of observations with a small error, and at the same time has small norm in the DFT domain.* This is done by solving a convex optimization problem (namely a well-structured SOCP).
The resulting estimator is shown to be adaptive to the unknown shift-invariant structure of the signal. Thus, the following signals are addressed:
- signals belonging to low-dimensional shift-invariant linear subspaces; equivalently, discretized solutions of low-order ODEs and PDEs, and in particular, sums of complex sinusoids;
- more generally, signals that are well-approximated in such subspaces; note that this includes the classical case of discretized smooth functions.

All the details regarding the theoretical motivation, statistical performance, and algorithmic implementation of the outlined approach can be found in [1,2,3]. We also plan to release the code reproducing the experimental results of those papers.

### Installation
Download or clone the repository, then simply change to the code directory.

### Demos
>> demo1d       % Run a demonstration of 1D-denoising
>> demo2d       % Run a demonstration of 2D-denoising

## Usage

## Documentation

### References

1. [Adaptive Recovery of Signals by Convex Optimization](https://hal.inria.fr/hal-01250215) Z. Harchaoui, A. Juditsky, A. Nemirovski, D. Ostrovskii
2. [Structure-Blind Signal Recovery](https://arxiv.org/abs/1607.05712) D. Ostrovskii, Z. Harchaoui, A. Judistky, A. Nemirovski
3. [Efficient First-Order Algorithms for Adaptive Signal Denoising](add a link) D. Ostrovskii, Z. Harchaoui
