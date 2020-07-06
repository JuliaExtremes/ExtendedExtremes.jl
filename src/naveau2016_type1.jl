using Distributions
import Distributions: cdf, sampler
import Random: AbstractRNG
import Base: rand

# macro for argument checking (récupéré de Distributions.jl)
macro check_args(D, cond)
    quote
        if !($(esc(cond)))
            throw(ArgumentError(string(
                $(string(D)), ": the condition ", $(string(cond)), " is not satisfied.")))
        end
    end
end

"""
    EGPpower(σ, ξ, κ)   (Nom provisoire pour la première famille)

Model (i), with the power law distribution G (v ) = vκ , is clearly the simplest choice and leads to a model in (6) with  three parameters: κ controls the shape of the lower tail, σ is a scale parameter, and ξ controls the rate of upper tail decay.

```julia
EGPpower(σ, ξ, κ)   # Extended Generalized Pareto of Naveau (type 1) with scale parameter σ, rate of upper tail decay ξ and shape of the lower tail κ.
```

Naveau, P., Huser, R., Ribereau, P., and Hannart, A. ( 2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, Water Resour. Res., 52, 2753– 2769, doi:10.1002/2015WR018552.
"""
struct EGPpower{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ::T    # rate of upper tail decay
    κ::T    # shape of the lower tail

    function EGPpower{T}(σ::T, ξ::T, κ::T) where {T <: Real}
        new{T}(σ, ξ, κ)
    end
end

function EGPpower(σ::T, ξ::T, κ::T; check_args=true) where {T <: Real}
    check_args && @check_args(EGPpower, σ > zero(σ) && κ > zero(κ))
    return EGPpower{T}(σ, ξ, κ)
end

EGPpower(σ::Real, ξ::Real, κ::Real) = EGPpower(promote(σ, ξ, κ)...)
EGPpower(σ::Integer, ξ::Integer, κ::Integer) = EGPpower(float(σ), float(ξ), float(κ))


minimum(d::EGPpower) = 0.0
maximum(d::EGPpower) = Inf
insupport(d::EGPpower, x::Real) = minimum(d) <= x <= maximum(d)


#### Parameters

scale(d::EGPpower) = d.σ
decay(d::EGPpower) = d.ξ    # noms à vérifier...
shape(d::EGPpower) = d.κ

params(d::EGPpower) = (d.σ, d.ξ, d.κ)
partype(::EGPpower{T}) where {T} = T


#### Evaluation

function logpdf(d::EGPpower{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, κ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    lg(v::Real) = log(κ) + (κ-1)*log(v)

    p = -log(σ) + Distributions.logpdf(pd, x/σ) + lg(Distributions.cdf(pd, x/σ))

    return p
end

pdf(d::EGPpower, x::Real) = exp(logpdf(d, x))

function logcdf(d::EGPpower{T}, x::Real) where T<:Real
    # À vérifier
    μ = 0
    (σ, ξ, κ) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    #G(v::Real, κ::Real) = v^κ
    #lG(v::Real, κ::Real) = κ*log(v)

    p = κ*Distributions.logcdf(pd, x/σ)

    return p
end

cdf(d::EGPpower, x::Real) = exp(logcdf(d, x))

function quantile(d::EGPpower{T}, p::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    (σ, ξ, κ) = params(d)

    invG(p::Real, κ::Real) = p^(1/κ)

    x = (σ/ξ)*(((1 - invG(p,κ))^(-ξ)) - 1)

    return x
end


#### Sampling

function rand(rng::AbstractRNG, d::EGPpower)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    invG(p::Real, κ::Real) = p^(1/κ)

    return (d.σ/d.ξ)*(((1 - invG(u,d.κ))^(-d.ξ)) - 1)
end

sampler(d::EGPpower) = d    # à vérifier
