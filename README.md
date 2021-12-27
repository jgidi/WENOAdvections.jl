# WENOAdvections.jl

In-development julia library to solve the advection equation in 1 dimension by means of the Weighted Essentially Non-Oscillatory (WENO) interpolation scheme
presented on Section 3 from [Level set equations on surfaces via de Closest Point Method (2008)](https://link.springer.com/article/10.1007/s10915-008-9196-6), by Colin Macdonald & Steven Ruuth.

Currently, only periodic boundary conditions are considered. Also, the solution of the advection equation in 2 dimensions is implemented by solving batches of 1-dimensional advections.

Proper documentation will be implemented in the near future. In the meantime, simple examples are provided within the `examples/` folder on the root of the repository.
