"""
    TNormal(κ)

"""
struct TNormal{T<:Real} <: ContinuousUnivariateDistribution
    κ::T
    TNormal{T}(κ::T) where {T<:Real} = new{T}(κ)
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