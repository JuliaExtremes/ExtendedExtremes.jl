# EGPD.jl*

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

**working title*, en attendant de trouver quelque chose pour extrêmes sous asymptotiques 😉


# Documentation de la librairie *EGPD.jl*

```julia
using EGPD
using DataFrames, Distributions, Gadfly;
```

## 1. Liste des features

### 1.1 Types de distributions

La librairie EGPD implémente les quatres familles paramétriques de Naveau et al. (2016) et la famille de Gamet et Jalbert (2021):

```julia
EGPpower,
EGPpowermix,
EGPbeta,
EGPbetapower,
EGPnormal
```

### 1.2 Méthodes

La librairie EGPD étend les méthodes de Distributions.jl. Voici les méthodes incluses dans la librairie :

```julia
cdf,         # cumulative distribution function
insupport,   # predicate, is x in the support of the distribution?
logpdf,      # log probability density
logcdf,      # cdf returning log-probability
maximum,
minimum,
params,      # get the tuple of parameters
partype,
pdf,         # probability density function
quantile,    # inverse of cdf (defined for p in (0,1))
rand,
sampler,
```

### 1.3 Estimation des paramètres

L'estimation des paramètres se fait par maximum de vraisemblance à l'aide des fonctions suivantes :

```julia
EGPpowerfit,
EGPpowermixfit,
EGPbetafit,
EGPbetapowerfit
```

## 2. Documentation

### 2.1 Types de distributions

#### 2.1.1 EGPpower


```julia
?EGPpower()
```


```
EGPpower(σ, ξ, κ)
```

*EGPpower* corresponds to the first extended GP model of Naveau et al. (2016), with the power law distribution G(v) = v^κ.

It is a three parameters family: κ controls the shape of the lower tail, σ is a scale parameter, and ξ controls the rate of upper tail decay.

```julia
EGPpower(σ, ξ, κ)   # EGP of Naveau et al. (2016) (type 1) with scale parameter σ, rate of upper tail decay ξ and shape of the lower tail κ.

params(d)           # Get the parameters, i.e. (σ, ξ, κ)
```

Reference :

  * Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.




#### 2.1.2 EGPpowermix


```julia
?EGPpowermix()
```




```
EGPpowermix(σ, ξ, κ₁, κ₂, p)
```

*EGPpowermix* corresponds to the second extended GP model of Naveau et al. (2016), which is a mixture of power laws with G(v) = pv^κ₁ + (1-p)v^κ₂.

It is a five parameters family: κ₁ controls the lower tail behavior, κ₂ modifies the shape of the density in its central part, σ is a scale parameter, and ξ is a shape parameter.

```julia
EGPpowermix(σ, ξ, κ₁, κ₂, p)   # EGP of Naveau et al. (2016) (type 2) with parameters σ, ξ, κ₁, κ₂, p.

params(d)           # Get the parameters, i.e. (σ, ξ, κ₁, κ₂, p)
```

Reference :

  * Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.




#### 2.1.3 EGPbeta


```julia
?EGPbeta()
```




```
EGPbeta(σ, ξ, δ)
```

*EGPbeta* corresponds to the third extended GP model of Naveau et al. (2016).

It is a three parameters family: δ is a threshold tuning parameter, σ is a scale parameter, and ξ is a shape parameter.

```julia
EGPbeta(σ, ξ, δ)   # EGP of Naveau et al. (2016) (type 3) with parameters σ, ξ, δ.

params(d)           # Get the parameters, i.e. (σ, ξ, δ)
```

Reference :

  * Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.




#### 2.1.4 EGPbetapower


```julia
?EGPbetapower()
```




```
EGPbetapower(σ, ξ, δ)
```

*EGPbetapower* corresponds to the fourth extended GP model of Naveau et al. (2016).

It is a four parameters family: δ is a threshold tuning parameter, κ is a parameter that controls the lower tail behavior, σ is a scale parameter, and ξ is a shape parameter.

```julia
EGPbetapower(σ, ξ, δ, κ)   # EGP of Naveau et al. (2016) (type 4) with parameters σ, ξ, δ, κ.

params(d)           # Get the parameters, i.e. (σ, ξ, δ, κ)
```

Reference :

  * Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.




### 2.2 Méthodes

#### 2.2.1 cdf


```julia

```

#### 2.2.2 insupport


```julia

```

#### 2.2.3 logpdf


```julia

```

#### 2.2.4 logcdf


```julia

```

#### 2.2.5 maximum


```julia

```

#### 2.2.6 minimum


```julia

```

#### 2.2.7 params


```julia

```

#### 2.2.8 partype


```julia

```

#### 2.2.9 pdf


```julia

```

#### 2.2.10 quantile


```julia

```

#### 2.2.11 rand


```julia

```

#### 2.2.12 sampler


```julia

```

### 2.3 Estimation des paramètres

#### 2.3.1 EGPpowerfit


```julia

```

#### 2.3.2 EGPpowermixfit


```julia

```

#### 2.3.3 EGPbetafit


```julia

```

#### 2.3.4 EGPbetapowerfit


```julia

```
