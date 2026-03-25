"""
    TBeta(α)

"""
struct TBeta{T<:Real, D} <: ContinuousUnivariateDistribution
    α::T
    _dist::D
    function TBeta{T}(α::T) where {T<:Real}
        a = one(T) / T(32)
        b = one(T) / T(2)
        dist = LocationScale(-a/(b-a), one(T)/(b-a), Truncated(Beta(α, α), a, b))
        return new{T, typeof(dist)}(α, dist)
    end
end

function TBeta(α::T; check_args=true) where {T <: Real}
    check_args && @check_args(TBeta, α > 0 )
    return TBeta{T}(α)
end

#### Outer constructors

TBeta() = TBeta(1.0, check_args=false)
TBeta(α::Int) = TBeta(float(α), check_args=false)

#### Parameters

params(pd::TBeta) = promote(pd.α)

#### Evaluations

minimum(::TBeta) = 0.0
maximum(::TBeta) = 1.0
insupport(pd::TBeta, x::Real) = minimum(pd) <= x <= maximum(pd)

getdistribution(pd::TBeta) = pd._dist

@inline cdf(pd::TBeta, x::Real) = cdf(pd._dist, x)
@inline logcdf(pd::TBeta, x::Real) = logcdf(pd._dist, x)
@inline logpdf(pd::TBeta, x::Real) = logpdf(pd._dist, x)
@inline pdf(pd::TBeta, x::Real) = pdf(pd._dist, x)
@inline quantile(pd::TBeta, p::Real) = quantile(pd._dist, p)
@inline rand(rng::AbstractRNG, pd::TBeta) = rand(rng, pd._dist)