# Release Notes

## 0.3.2

- Remove the use of the `Distributions.@check_args` macro from constructors to ensure compatibility with `Distributions.jl` versions later than v0.25.126.

## 0.3.1

- Corrected the censored likelihood function.
- Updated the EGP likelihood to use log κ and log σ.
- Add an optional left censoring bound for the fit of the extended genralized Pareto distribution.
