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

function getdistribution(pd::Power)
   
    κ = params(pd)[1]
    
    return Beta(κ, 1)
    
end

function cdf(pd::Power, x::Real)
   
    td = getdistribution(pd)
    
    return cdf(td, x)
    
end

function logpdf(pd::Power, x::Real)
   
    td = getdistribution(pd)
    
    return logpdf(td, x)
    
end

function pdf(pd::Power, x::Real)
   
    td = getdistribution(pd)
    
    return pdf(td, x)
    
end

function quantile(pd::Power, p::Real)
    
    td = getdistribution(pd)
    
    return quantile(td, p)
    
end

function rand(rng::AbstractRNG, pd::Power)
    
    td = getdistribution(pd)

    return rand(rng, td)
end