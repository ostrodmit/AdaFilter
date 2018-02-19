# AdaFilter — Adaptive Denoising via Non-Linear Filtering

AdaFilter is a **MATLAB** package which implements **struture-adaptive denoising** of discrete-time signals and images.

## Description

The approach is to fit from the observations a time-invariant filter which *reproduces the vector of observations with a small error, and at the same time has small norm in the DFT domain.* This is done by solving a certain convex optimization problem, see [3].
Such estimators are adaptive to the unknown shift-invariant structure of the signal. Thus, the following cases are addressed:

- signals belonging to low-dimensional shift-invariant linear subspaces; equivalently, discretized solutions of low-order ODEs and PDEs, and in particular, sums of complex sinusoids;
- more generally, signals that are well-approximated in such subspaces; this includes the classical case of discretized smooth functions.

All the details regarding the theoretical motivation, statistical performance, and algorithmic implementation of the outlined approach can be found in [1,2,3]. We also plan to release the code reproducing the experimental results of those papers.

## Installation
Download or clone the repository, then change to the code directory: 
```
>> cd adafilter-master/code
``` 

## Usage
The user deals with two MATLAB functions, ``filter_recovery`` and ``fit_filter``, implemented in eponymous .m-files.

- ``filter_recovery`` is a high-level interface for denoising. It accepts three input parameters: 
  - 1D or 2D array of observations ``y``, a noisy version of the signal ``x``
  - structure ``params`` with denoising parameters; 
  - optionally, structure ``solver_control`` with parameters of the filter-fitting procedure which are passed to ``fit_filter``.
The output is an estimate of the input signal. 

- ``fit_filter`` contains implementations of the filter-fitting procedures — specific optimization problems, and first-order solvers for them, see [3].
*This function must only be accessed if fine-tuning of the approach is required.*
It accepts two arrays ``y1`` and ``y2``, and computes a one-sided filter which reproduces ``y2`` via the convolution of ``y1`` with the filter. The third arguument is the structure ``control`` with the parameters of a filter-fitting procedure (see the built-in documentation).

## Documentation
The documentation of both functions is available via the ``doc`` command. In **MATLAB**, run
```
>> doc <function>
```
where ``<function>`` is ``filter_recovery`` or ``fit_filter``. The following conventions are used: 
- The obligatory parameters are marked by ``*``.
- A parameter is sometimes succeeded by ``[<Range>,<Default>]``. In this case, ``<Range>``is a symbolic description of the set of *admissible values* of a parameter, and ``<Default>`` specifies the *default value* assigned if the corresponding field of the input structure is empty. f it is replaced with ``*``, a description in words follows later.

## Demos
We provide "numerical tours" to demonstrate application of our approach to different signals in 1D (``demo1d.m``) and 2D (``demo2d.m``). E.g. for the 1D demo, run
```
publish('demo1d.m',<format>)
```
where ``<format>`` is the format in which the manual file will be generated; good choices are ``'.html'`` and ``'.pdf'``.

## Features
As of Feb. 2018, the following features are implemented:
- Lepski-type bandwidth adaptation;
- pointwise and blockwise filtering mode;
- parallelized computation via ``parfor``.

For their description, see the documentation of ``filter_recovery``.

## References
1. [Adaptive Recovery of Signals by Convex Optimization](https://hal.inria.fr/hal-01250215) Z. Harchaoui, A. Juditsky, A. Nemirovski, D. Ostrovskii
2. [Structure-Blind Signal Recovery](https://arxiv.org/abs/1607.05712) D. Ostrovskii, Z. Harchaoui, A. Judistky, A. Nemirovski
3. [Efficient First-Order Algorithms for Adaptive Signal Denoising](https://arxiv.org/abs/1607.05712) D. Ostrovskii, Z. Harchaoui
