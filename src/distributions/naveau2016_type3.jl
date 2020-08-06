"""
    EGPbeta(σ, ξ, δ)

*EGPbeta* corresponds to the third extended GP model of Naveau et al. (2016).

It is a three parameters family: δ is a threshold tuning parameter, σ is a scale parameter, and ξ is a shape parameter.

```julia
EGPbeta(σ, ξ, δ)   # EGP of Naveau et al. (2016) (type 3) with parameters σ, ξ, δ.

params(d)           # Get the parameters, i.e. (σ, ξ, δ)
```

Reference :

* Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.
"""
struct EGPbeta{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ::T    # shape parameter
    δ::T    # threshold tuning parameter

    function EGPbeta{T}(σ::T, ξ::T, δ::T) where {T <: Real}
        new{T}(σ, ξ, δ)
    end
end

function EGPbeta(σ::T, ξ::T, δ::T; check_args=true) where {T <: Real}
    check_args && @check_args(EGPbeta, σ > zero(σ) && δ > zero(δ))
    return EGPbeta{T}(σ, ξ, δ)
end

EGPbeta(σ::Real, ξ::Real, δ::Real) = EGPbeta(promote(σ, ξ, δ)...)
EGPbeta(σ::Integer, ξ::Integer, δ::Integer) = EGPbeta(float(σ), float(ξ), float(δ))


minimum(d::EGPbeta) = 0.0
maximum(d::EGPbeta) = Inf * (d.ξ >= 0) - ( d.σ / d.ξ ) * (d.ξ < 0)
insupport(d::EGPbeta, x::Real) = minimum(d) <= x <= maximum(d)


#### Parameters

scale(d::EGPbeta) = d.σ
shape(d::EGPbeta) = d.ξ    # noms à vérifier...
thresh(d::EGPbeta) = d.δ

params(d::EGPbeta) = (d.σ, d.ξ, d.δ)
partype(::EGPbeta{T}) where {T} = T


#### Evaluation

function logpdf(d::EGPbeta{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, δ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    #g(v::Real, δ) = δ*(1-v)^(δ-1) * pdf(Beta(1/δ, 2), (1 - v)^δ)
    lg(v::Real, δ) = log(δ) + (δ-1)*log(1-v) + Distributions.logpdf(Beta(1/δ, 2), (1 - v)^δ)

    lf = -log(σ) + Distributions.logpdf.(pd, x/σ) + lg(Distributions.cdf(pd, x/σ), δ)

    return lf
end

pdf(d::EGPbeta, x::Real) = exp(logpdf(d, x))

function logcdf(d::EGPbeta{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, δ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    G(v::Real, δ::Real) = 1 - Distributions.cdf(Beta(1/δ, 2), (1 - v)^δ)

    lf = log(G(Distributions.cdf(pd, x/σ), δ))

    return lf
end

cdf(d::EGPbeta, x::Real) = exp(logpdf(d, x))

function quantile(d::EGPbeta{T}, p::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    (σ, ξ, δ) = params(d)

    invG(v::Real, δ::Real) = 1 - Distributions.quantile(Beta(1/δ, 2), (1 - v))^(1/δ)

    x = (σ/ξ)*(((1 - invG(p,δ))^(-ξ)) - 1)

    return x
end


#### Sampling

function rand(rng::AbstractRNG, d::EGPbeta)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    invG(v::Real, δ::Real) = 1 - Distributions.quantile(Beta(1/δ, 2), (1 - v))^(1/δ)

    return (d.σ/d.ξ)*(((1 - invG(u,d.δ))^(-d.ξ)) - 1)
end

sampler(d::EGPbeta) = d
