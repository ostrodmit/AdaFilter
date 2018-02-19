# AdaFilter — Signal Recovery via Adaptive Filtering

AdaFilter provides a MATLAB implementation of struture-adaptive signal and image denoising.

## Short description

The goal is to estimate a discrete-time signal 

In a nutshell, one fits a data-dependent time-invariant filter by solving an appropriate (well-structured) convex optimization problem. 
This approach allows to adapt to the unknown shift-invariant structure of the signal, encompassing, in particular, signals satisfying low-order ODEs and PDEs, and more generally, signals 

All the details concerning the theore
