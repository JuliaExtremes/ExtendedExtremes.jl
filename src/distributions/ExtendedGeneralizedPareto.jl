"""
    ExtendedGeneralizedPareto()

"""
struct ExtendedGeneralizedPareto{T} <: ContinuousUnivariateDistribution
    # Distribution on (0,1)
    V::T
    # Tail distribution
    G::GeneralizedPareto
end

ExtendedGeneralizedPareto() = ExtendedGeneralizedPareto(Uniform(), GeneralizedPareto())

#TODO: Verify V is in (0,1)

#### Parameters

scale(pd::ExtendedGeneralizedPareto) = scale(pd.G)
tailindex(pd::ExtendedGeneralizedPareto) = shape(pd.G)
shape(pd::ExtendedGeneralizedPareto) = params(pd.V)

params(pd::ExtendedGeneralizedPareto) = (params(pd.V)..., params(pd.G)...)
#partype(pd::ExtendedGeneralizedPareto) where {T} = T

#### Evaluations

function EGPtype(pd::Type{<:ExtendedGeneralizedPareto{T}}) where T<:ContinuousUnivariateDistribution
    return T
end

minimum(pd::ExtendedGeneralizedPareto) = minimum(pd.G)
maximum(pd::ExtendedGeneralizedPareto) = maximum(pd.G)
insupport(pd::ExtendedGeneralizedPareto, x::Real) = minimum(pd) <= x <= maximum(pd)

function cdf(pd::ExtendedGeneralizedPareto, x::Real)
   
    b = cdf(pd.G, x)
    
    return cdf(pd.V, b)
    
end

function logpdf(pd::ExtendedGeneralizedPareto, x::Real)
    
    b = cdf(pd.G, x)

    logdensity = logpdf(pd.V, b) + logpdf(pd.G, x)
       
    return logdensity
    
end

function pdf(pd::ExtendedGeneralizedPareto, x::Real)
   
    return exp(logpdf(pd, x))
    
end

function quantile(pd::ExtendedGeneralizedPareto, p::Real)
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    b = quantile(pd.V, p)
    
    return quantile(pd.G, b)

end