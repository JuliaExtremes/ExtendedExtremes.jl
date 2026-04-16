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
    x <= 0. && return 0.0
    x >= 1. && return 1.0
    return x^pd.κ
end

function logcdf(pd::Power, x::Real)
    x <= 0. && return -Inf
    x >= 1 && return 0.0
    return pd.κ * log(x)
end

function logpdf(pd::Power, x::Real)
    κ = pd.κ
    (0. < x < 1.) && return -Inf
    x == 0. && return κ > 1. ? -Inf : (κ < 1. ? Inf : 0.0)
    return log(κ) + (κ - 1) * log(x)
end

function pdf(pd::Power, x::Real)
    κ = pd.κ
    (0. < x < 1.) && return 0.0
    x == 0. && κ > 1. && return 0.0
    return κ * x^(κ - 1.)
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
    return rand(rng)^(1. / pd.κ)
end