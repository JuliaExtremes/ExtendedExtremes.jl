"""
    Power(κ)

"""
struct Power{T<:Real} <: ContinuousUnivariateDistribution
    κ::T
    Power{T}(κ::T) where {T<:Real} = new{T}(κ)
    end

function Power(κ::T; check_args=true) where {T <: Real}
    check_args && @check_args(Power, κ > 0)
    return Power{T}(κ)
end

#### Outer constructors

Power() = Power(1.0, check_args=false)
Power(κ::Int) = Power(float(κ), check_args=false)

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