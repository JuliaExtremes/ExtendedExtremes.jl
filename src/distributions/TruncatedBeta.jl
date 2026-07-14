"""
    TBeta(α)

"""
struct TBeta{T<:Real} <: ContinuousUnivariateDistribution
    α::T

    function TBeta{T}(α::T) where {T<:Real}
        α > zero(α) || throw(DomainError(α, "TBeta: α must be strictly positive."))

        return new{T}(α)
    end
end

#### Outer constructors

TBeta() = TBeta(1.0)
TBeta(α::Int) = TBeta(float(α))

#### Parameters

params(pd::TBeta) = promote(pd.α)

#### Evaluations

minimum(::TBeta) = 0.0
maximum(::TBeta) = 1.0
insupport(pd::TBeta, x::Real) = minimum(pd) <= x <= maximum(pd)

function getdistribution(pd::TBeta)
   
    α = params(pd)[1]
    
    a = 1/32
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