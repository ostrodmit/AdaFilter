# AdaFilter — Signal Denoising via Adaptive Filtering

AdaFilter provides a **MATLAB** implementation of struture-adaptive denoising of discrete-time signals and images.

## Description

The approach is to fit from the observations a time-invariant filter which *reproduces the vector of observations with a small error, and at the same time has small norm in the DFT domain.* This is done by solving a convex optimization problem (namely a well-structured SOCP).
The resulting estimator is shown to be adaptive to the unknown shift-invariant structure of the signal. Thus, the following cases are addressed:

- signals belonging to low-dimensional shift-invariant linear subspaces; equivalently, discretized solutions of low-order ODEs and PDEs, and in particular, sums of complex sinusoids;
- more generally, signals that are well-approximated in such subspaces; this includes the classical case of discretized smooth functions.

All the details regarding the theoretical motivation, statistical performance, and algorithmic implementation of the outlined approach can be found in [1,2,3]. We also plan to release the code reproducing the experimental results of those papers.

## Installation
Download or clone the repository, then change to the code directory: 
```
>> cd adafilter-master/code
``` 

## Usage and Documentation
The user deals with two MATLAB functions, ``filter_recovery`` and ``fit_filter``, implemented in eponymous .m-files.

- ``filter_recovery`` is a high-level interface for denoising. It accepts three input parameters: 
  - 1D or 2D array of observations ``y`` — a noisy version of the signal ``x``
  - structure ``params`` with denoising parameters; 
  - optionally, structure ``solver_control`` with parameters of the filter-fitting procedure which are passed to ``fit_filter``.
The output is an estimate of the input signal. To see the documentation and example of usage, run 
``
>> doc filter_recovery
``

- ``fit_filter`` 

The documentation of both functions is available from inside these functions via the ``doc`` command.

## Features
As of Feb. 2018, the following features are implemented:

## Demos
We provide "numerical tours" to demonstrate application of our approach to different signals in 1D and 2D. Follow these steps:
1. Run **MATLAB** GUI with administrative rights. 
2. Open ``demo1d.m`` or ``demo2d.m``.
3. Go to **Publish** tab in the main menu, and press **Publish**. The options can be edited in the dropdown menu.
MATLAB will generate an .html-file, and automatically open it in the default browser.

## References
1. [Adaptive Recovery of Signals by Convex Optimization](https://hal.inria.fr/hal-01250215) Z. Harchaoui, A. Juditsky, A. Nemirovski, D. Ostrovskii
2. [Structure-Blind Signal Recovery](https://arxiv.org/abs/1607.05712) D. Ostrovskii, Z. Harchaoui, A. Judistky, A. Nemirovski
3. [Efficient First-Order Algorithms for Adaptive Signal Denoising](add a link) D. Ostrovskii, Z. Harchaoui
