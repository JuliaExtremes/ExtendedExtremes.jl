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