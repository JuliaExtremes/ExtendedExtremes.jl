"""
    Power(κ)

"""
struct Power{T<:Real} <: ContinuousUnivariateDistribution
    κ::T

    function Power{T}(κ::T) where {T<:Real}
        κ > zero(κ) || throw(DomainError(κ, "Power: κ must be strictly positive."))
        return new{T}(κ)
    end
end

#### Outer constructors

Power(κ::T) where {T<:AbstractFloat} = Power{T}(κ)
Power(κ::Integer) = Power(float(κ))
Power() = Power(1.0)

#### Parameters

params(pd::Power) = promote(pd.κ)

#### Evaluations

minimum(::Power) = 0.0
maximum(::Power) = 1.0
insupport(pd::Power, x::Real) = minimum(pd) <= x <= maximum(pd)

function cdf(pd::Power, x::Real)
    return exp(logcdf(pd, x))
end

function logcdf(pd::Power, x::Real)
    temp = x < zero(x) ? oftype(x, -Inf) : pd.κ*log(x)
    return x > one(x) ? zero(x) : temp
end

function logpdf(pd::Power, x::Real)
    κ = pd.κ
    p = log(κ) + (κ - 1.)*log(x)
    return (zero(x) < x < one(x)) ? p : oftype(p, -Inf)
end

function loglikelihood(pd::Power, x::Vector{<:Real})
    κ = pd.κ
    n = length(x)
    s = sum(log, x)
    return n*log(κ) + (κ - 1.)*s 
end

function quantile(pd::Power, p::Real)
    @assert zero(p) < p < one(p)
    return p^(1. / pd.κ)    
end

function rand(rng::AbstractRNG, pd::Power)
    dist = Beta(pd.κ, 1.)
    return rand(rng, dist)
end