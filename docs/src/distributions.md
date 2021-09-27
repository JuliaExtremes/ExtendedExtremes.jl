# Extended GP Distributions

## Distributions

The extended generalized Pareto families proposed in Naveau et al. (2016) and the extended GP model developed by Gamet and Jalbert (2020) are implemented and ready to be used just like any [`UnivariateDistribution`](https://juliastats.org/Distributions.jl/stable/univariate/) from the *[Distributions.jl](https://github.com/JuliaStats/Distributions.jl)* package.

```@docs
EGPpower
EGPpowermix
EGPbeta
EGPbetapower
EGPnormal
```

## Common Interface

A series of methods are implemented for each EGPD, which provide useful functionalities such as moment computation, pdf evaluation, and sampling (i.e. random number generation).

### Probability Evaluation

```@docs
cdf(d::UnivariateDistribution, x::Real)
logcdf(d::UnivariateDistribution, x::Real)
logpdf(d::UnivariateDistribution, x::Real)
pdf(d::UnivariateDistribution, x::Real)
quantile(d::UnivariateDistribution, q::Real)
```

### Computation of statistics

```@docs
insupport(d::UnivariateDistribution, x::Any)
maximum(d::UnivariateDistribution)
minimum(d::UnivariateDistribution)
```

### Parameter Retrieval

```@docs
params(d::UnivariateDistribution)
```

## Index

```@index
Pages = ["distributions.md"]
```