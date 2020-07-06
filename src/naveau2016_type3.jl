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
maximum(d::EGPbeta) = Inf
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
