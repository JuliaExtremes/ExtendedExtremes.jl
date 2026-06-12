# Release Notes

## Nightly

- Improve computation performance of `logpdf()` for the TBeta distribution. It is not relying anymore on `Truncated()` and `LocationScale()` by implementing directly the log-density.   

## 0.3.2

- Corrected the censored likelihood function.
- Updated the EGP likelihood to use log κ and log σ.

## 0.3.1

- Add an optional left censoring bound for the fit of the extended genralized Pareto distribution.