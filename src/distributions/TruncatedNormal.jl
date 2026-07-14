"""
    TNormal(κ)

"""
struct TNormal{T<:Real} <: ContinuousUnivariateDistribution
    κ::T

    function TNormal{T}(κ::T) where {T<:Real}
        κ > zero(κ) || throw(DomainError(κ, "TNormal: κ must be strictly positive."))

        return new{T}(κ)
    end
end

#### Outer constructors

TNormal(κ::T) where {T<:Real} = TNormal{T}(κ)
TNormal() = TNormal(1.0)
TNormal(κ::Int) = TNormal(float(κ))

#### Parameters

params(pd::TNormal) = promote(pd.κ)

#### Evaluations

minimum(::TNormal) = 0.0
maximum(::TNormal) = 1.0
insupport(pd::TNormal, x::Real) = minimum(pd) <= x <= maximum(pd)

function getdistribution(pd::TNormal)
   
    κ = params(pd)[1]
    
    return Truncated(Normal(1, sqrt(1/κ)), 0, 1)
    
end

function cdf(pd::TNormal, x::Real)
   
    td = getdistribution(pd)
    
    return cdf(td, x)
    
end

function logpdf(pd::TNormal, x::Real)
   
    td = getdistribution(pd)
    
    return logpdf(td, x)
    
end

function pdf(pd::TNormal, x::Real)
   
    td = getdistribution(pd)
    
    return pdf(td, x)
    
end

function quantile(pd::TNormal, p::Real)
    
    td = getdistribution(pd)
    
    return quantile(td, p)
    
end

function rand(rng::AbstractRNG, pd::TNormal)
    
    td = getdistribution(pd)

    return rand(rng, td)
end