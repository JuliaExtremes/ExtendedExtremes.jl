"""
    TNormal(κ)

"""
struct TNormal{T<:Real, D} <: ContinuousUnivariateDistribution
    κ::T
    _dist::D
    function TNormal{T}(κ::T) where {T<:Real}
        dist = Truncated(Normal(one(T), sqrt(one(T)/κ)), zero(T), one(T))
        return new{T, typeof(dist)}(κ, dist)
    end
end

function TNormal(κ::T; check_args=true) where {T <: Real}
    check_args && @check_args(TNormal, κ > 0)
    return TNormal{T}(κ)
end

#### Outer constructors

TNormal() = TNormal(1.0, check_args=false)
TNormal(κ::Int) = TNormal(float(κ), check_args=false)

#### Parameters

params(pd::TNormal) = promote(pd.κ)

#### Evaluations

minimum(::TNormal) = 0.0
maximum(::TNormal) = 1.0
insupport(pd::TNormal, x::Real) = minimum(pd) <= x <= maximum(pd)

getdistribution(pd::TNormal) = pd._dist

@inline cdf(pd::TNormal, x::Real) = cdf(pd._dist, x)
@inline logpdf(pd::TNormal, x::Real) = logpdf(pd._dist, x)
@inline pdf(pd::TNormal, x::Real) = pdf(pd._dist, x)
@inline quantile(pd::TNormal, p::Real) = quantile(pd._dist, p)
@inline rand(rng::AbstractRNG, pd::TNormal) = rand(rng, pd._dist)