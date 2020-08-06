"""
    EGPbetapower(σ, ξ, δ)

*EGPbetapower* corresponds to the fourth extended GP model of Naveau et al. (2016).

It is a four parameters family: δ is a threshold tuning parameter, κ is a parameter that controls the lower tail behavior, σ is a scale parameter, and ξ is a shape parameter.

```julia
EGPbetapower(σ, ξ, δ, κ)   # EGP of Naveau et al. (2016) (type 4) with parameters σ, ξ, δ, κ.

params(d)           # Get the parameters, i.e. (σ, ξ, δ, κ)
```

Reference :

* Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.
"""
struct EGPbetapower{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ::T    # shape parameter
    δ::T    # threshold tuning parameter
    κ::T    # lower tail parameter

    function EGPbetapower{T}(σ::T, ξ::T, δ::T, κ::T) where {T <: Real}
        new{T}(σ, ξ, δ, κ)
    end
end

function EGPbetapower(σ::T, ξ::T, δ::T, κ::T; check_args=true) where {T <: Real}
    check_args && @check_args(EGPbetapower, σ > zero(σ) && δ > zero(δ) && κ > zero(κ))
    return EGPbetapower{T}(σ, ξ, δ, κ)
end

EGPbetapower(σ::Real, ξ::Real, δ::Real, κ::Real) = EGPbetapower(promote(σ, ξ, δ, κ)...)
EGPbetapower(σ::Integer, ξ::Integer, δ::Integer, κ::Integer) = EGPbetapower(float(σ), float(ξ), float(δ), float(κ))


minimum(d::EGPbetapower) = 0.0
maximum(d::EGPbetapower) = Inf * (d.ξ >= 0) - ( d.σ / d.ξ ) * (d.ξ < 0)
insupport(d::EGPbetapower, x::Real) = minimum(d) <= x <= maximum(d)


#### Parameters

scale(d::EGPbetapower) = d.σ
shape(d::EGPbetapower) = d.ξ    # noms à vérifier...
thresh(d::EGPbetapower) = d.δ
lowertail(d::EGPbetapower) = d.κ

params(d::EGPbetapower) = (d.σ, d.ξ, d.δ, d.κ)
partype(::EGPbetapower{T}) where {T} = T


#### Evaluation

function logpdf(d::EGPbetapower{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, δ, κ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)
    beta_pd = Beta(1/δ, 2)

    #g(v::Real, δ::Real, κ::Real) = (κ/2)*(1 - cdf(Beta(1/δ,2), (1-v)^δ))^(κ/2 - 1) * δ*(1-v)^(δ-1) * pdf(Beta(1/δ, 2), (1 - v)^δ)
    lg(v::Real, δ::Real, κ::Real) = log(κ/2) + (κ/2 - 1)*log(1 - Distributions.cdf(Beta(1/δ,2), (1-v)^δ)) + log(δ) + (δ-1)*log(1-v) + Distributions.logpdf(beta_pd, (1 - v)^δ)

    lf = -log(σ) + Distributions.logpdf.(pd, x/σ) + lg(Distributions.cdf(pd, x/σ), δ, κ)  # log de la vraisemblance

    return lf
end

pdf(d::EGPbetapower, x::Real) = exp(logpdf(d, x))

function logcdf(d::EGPbetapower{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, δ, κ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    G(v::Real, δ::Real, κ::Real) = (1 - Distributions.cdf(Beta(1/δ, 2), (1 - v)^δ))^(κ/2)

    lf = log(G(Distributions.cdf(pd, x/σ), δ, κ))

    return lf
end

cdf(d::EGPbetapower, x::Real) = exp(logpdf(d, x))

function quantile(d::EGPbetapower{T}, p::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    (σ, ξ, δ, κ) = params(d)

    invG(v::Real, δ::Real, κ::Real) = 1 - Distributions.quantile(Beta(1/δ, 2), (1 - v^(2/κ)))^(1/δ)

    x = (σ/ξ)*(((1 - invG(p,δ,κ))^(-ξ)) - 1)

    return x
end


#### Sampling

function rand(rng::AbstractRNG, d::EGPbetapower)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    invG(v::Real, δ::Real, κ::Real) = 1 - Distributions.quantile(Beta(1/δ, 2), (1 - v^(2/κ)))^(1/δ)

    return (d.σ/d.ξ)*(((1 - invG(u,d.δ,d.κ))^(-d.ξ)) - 1)
end

sampler(d::EGPbetapower) = d
