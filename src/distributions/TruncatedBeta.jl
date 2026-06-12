"""
    TBeta(α)

"""
struct TBeta{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    TBeta{T}(α::T) where {T<:Real} = new{T}(α)
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


function cdf(pd::TBeta, x::Real)
    κ = pd.α

    if x <= 0
        return 0.0
    elseif x >= 1
        return 1.0
    end

    a = 1 / 32
    b = 1 / 2
    width = b - a

    w = a + width * x

    dist = Beta(κ, κ)

    Fa = cdf(dist, a)
    Fw = cdf(dist, w)

    Z = 1 / 2 - Fa

    return (Fw - Fa) / Z
end

function logpdf(pd::TBeta, x::Real)
    κ = pd.α

    if x < 0 || x > 1
        return -Inf
    end

    a = 1 / 32
    width = 1/2 - a

    w = a + width * x

    dist = Beta(κ, κ)

    Fa = cdf(dist, a)
    Z = 1 / 2 - Fa

    return log(width) + logpdf(dist , w) - log(Z)
end

function quantile(pd::TBeta, p::Real)
    κ = pd.α

    if p < 0 || p > 1
        throw(ArgumentError("p must be in [0, 1]"))
    elseif p == 0
        return 0.0
    elseif p == 1
        return 1.0
    end

    a = 1 / 32
    width = 1/2 - a

    dist = Beta(κ, κ)

    Fa = cdf(dist, a)
    Fb = 1 / 2

    w = quantile(dist, Fa + p * (Fb - Fa))

    return (w - a) / width
end


function rand(rng::AbstractRNG, pd::TBeta)
    u = rand(rng)
    return quantile(pd, u)
end

rand(pd::TBeta) = rand(Random.default_rng(), pd)