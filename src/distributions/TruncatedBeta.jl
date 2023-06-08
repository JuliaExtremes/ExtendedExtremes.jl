"""
    TBeta(α)

"""
struct TBeta{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    α::T
    TBeta{T}(a::T, α::T) where {T<:Real} = new{T}(a,α)
    end

function TBeta(a::T, α::T; check_args=true) where {T <: Real}
    check_args && @check_args(TBeta, 0<a<1/2,  α > 0 )
    return TBeta{T}(a, α)
end

#### Outer constructors

TBeta(α::Real=1.0) = TBeta(1/32, float(α), check_args=true)
TBeta(α::Int) = TBeta(1/32, float(α), check_args=true)
TBeta(a::Real, α::Real) = TBeta(promote(a, α)..., check_args=true)


#### Parameters

params(pd::TBeta) = (pd.a, pd.α)

#### Evaluations

minimum(::TBeta) = 0.0
maximum(::TBeta) = 1.0
insupport(pd::TBeta, x::Real) = minimum(pd) <= x <= maximum(pd)

function getdistribution(pd::TBeta)
   
    a, α = params(pd)

    b = 1/2
    
    return LocationScale(-a/(b-a), 1/(b-a), Truncated(Beta(α, α), a, b))
    
end

function cdf(pd::TBeta, x::Real)
   
    td = getdistribution(pd)
    
    return cdf(td, x)
    
end

function logpdf(pd::TBeta, x::Real)
   
    td = getdistribution(pd)
    
    return logpdf(td, x)
    
end

function pdf(pd::TBeta, x::Real)
   
    td = getdistribution(pd)
    
    return pdf(td, x)
    
end

function quantile(pd::TBeta, p::Real)
    
    td = getdistribution(pd)
    
    return quantile(td, p)
    
end

function rand(rng::AbstractRNG, pd::TBeta)
    
    td = getdistribution(pd)

    return rand(rng, td)
end