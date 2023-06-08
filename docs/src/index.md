# ExtendedExtremes Package

*[ExtendedExtremes.jl](https://github.com/JuliaExtremes/ExtendedExtremes.jl)* provides a collection of extended Generalized Pareto distributions (EGPD) and related functions. In particular, the extended generalized Pareto distributions described in [Gamet & Jalbert (2022)](https://doi.org/10.1002/env.2744) and related methods are implemented, such as:

* Sampling from distributions
* Probability density functions (pdf) and their logarithm (logpdf)
* Cumulative distribution functions (cdf) and their logarithm (logcdf)
* Maximum likelihood estimation
* Diagnostic plots for assessing model accuracy

!!! note

	The *ExtendedExtremes.jl* library is currently being refactored. The 	support for extended generalized Pareto families proposed in Naveau et al. 	(2016) has been abandoned as well as some other functionalities. However, it is still possible to find these abandoned features in *ExtendedExtremes.jl* 	[v0.1.1](https://github.com/JuliaExtremes/ExtendedExtremes.jl/releases/tag/	v0.1.1).



### Reference

Gamet, P. & Jalbert, J. (2022). A flexible extended generalized Pareto distribution for tail estimation. *Environmetrics*, 33(6), e2744.