"""
  EGPnormal(σ, ξ, κ)
*EGPnormal* corresponds to the extended GP model developed by Gamet and Jalbert (2020), with the truncated normal distribution of mean 1 and variance 1/κ^2.
It is a three parameters family: κ controls the shape of the lower tail, σ is a scale parameter, and ξ controls the rate of upper tail decay.
```julia
EGPnormal(σ, ξ, κ)   # EGP of Gamet and Jalbert (2020) with scale parameter σ, rate of upper tail decay ξ and shape of the lower tail κ.
params(d)           # Get the parameters, i.e. (σ, ξ, κ)
```
Reference :
unpublished for now
"""
struct EGPnormal{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ::T    # rate of upper tail decay
    κ::T    # shape of the lower tail

    function EGPnormal{T}(σ::T, ξ::T, κ::T) where {T <: Real}
        new{T}(σ, ξ, κ)
    end
end

function EGPnormal(σ::T, ξ::T, κ::T; check_args=true) where {T <: Real}
    check_args && @check_args(EGPnormal, σ > zero(σ) && κ >= zero(κ))
    return EGPnormal{T}(σ, ξ, κ)
end

EGPnormal(σ::Real, ξ::Real, κ::Real) = EGPnormal(promote(σ, ξ, κ)...)
EGPnormal(σ::Integer, ξ::Integer, κ::Integer) = EGPnormal(float(σ), float(ξ), float(κ))


minimum(d::EGPnormal) = 0.0
maximum(d::EGPnormal) = Inf * d.ξ >= 0 - ( d.σ / d.ξ ) * d.ξ < 0
insupport(d::EGPnormal, x::Real) = minimum(d) <= x <= maximum(d)


#### Parameters

scale(d::EGPnormal) = d.σ
decay(d::EGPnormal) = d.ξ    # noms à vérifier...
shape(d::EGPnormal) = d.κ

params(d::EGPnormal) = (d.σ, d.ξ, d.κ)
partype(::EGPnormal{T}) where {T} = T


#### Evaluation

function pdf(d::EGPnormal{T}, x::Real) where T<:Real

    (σ, ξ, κ) = params(d)
    η = 1

    if κ == 0

        return Distributions.pdf.(GeneralizedPareto(0, σ, ξ), x)

    else

        return κ/(Distributions.cdf(Normal(), κ*(1-η))-Distributions.cdf(Normal(), -κ*η))*Distributions.pdf(Normal(), κ*(Distributions.cdf(GeneralizedPareto(0, σ, ξ),x)-η))*Distributions.pdf(GeneralizedPareto(0, σ, ξ),x)

    end
end

logpdf(d::EGPnormal, x::Real) = log(pdf(d, x))

function cdf(d::EGPnormal{T}, x::Real) where T<:Real

    (σ, ξ, κ) = params(d)

    return 2/Distributions.erf(κ/sqrt(2)) * (Distributions.cdf(Normal(), κ*(Distributions.cdf(GeneralizedPareto(σ, ξ), x)-1)) - 1 + Distributions.cdf(Normal(), κ))

end

logcdf(d::EGPnormal, x::Real) = log(cdf(d, x))

function quantile(d::EGPnormal{T}, p::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    (σ, ξ, κ) = params(d)

    return Distributions.quantile(GeneralizedPareto(σ,ξ), 1 + 1/κ * Distributions.quantile(Normal(), 0.5*(1 - (1-p)*Distributions.erf(κ/sqrt(2)))))

end


#### Sampling

function rand(rng::AbstractRNG, d::EGPnormal)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    if u == 1

        return maximum(d)

    else

        return quantile(d, u)

    end
end

sampler(d::EGPnormal) = d
