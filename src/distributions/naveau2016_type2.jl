struct EGPpowermix{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ::T    # rate of upper tail decay
    κ₁::T    # shape of the lower tail
    κ₂::T   # shape of the density in the central part
    p::T    # probability

    function EGPpowermix{T}(σ::T, ξ::T, κ₁::T,κ₂::T, p::T) where {T <: Real}
        new{T}(σ, ξ, κ₁, κ₂, p)
    end
end

function EGPpowermix(σ::T, ξ::T, κ₁::T, κ₂::T, p::T; check_args=true) where {T <: Real}
    check_args && @check_args(EGPpowermix, σ > zero(σ) && κ₁ > zero(κ₁) && κ₂ > zero(κ₂))  # manque test pour p
    return EGPpowermix{T}(σ, ξ, κ₁, κ₂, p)
end

EGPpowermix(σ::Real, ξ::Real, κ₁::Real, κ₂::Real, p::Real) = EGPpowermix(promote(σ, ξ, κ₁, κ₂, p)...)
EGPpowermix(σ::Integer, ξ::Integer, κ₁::Integer, κ₂::Integer, p::Integer) = EGPpowermix(float(σ), float(ξ), float(κ₁), float(κ₂), float(p))


minimum(d::EGPpowermix) = 0.0
maximum(d::EGPpowermix) = Inf
insupport(d::EGPpowermix, x::Real) = minimum(d) <= x <= maximum(d)


#### Parameters

scale(d::EGPpowermix) = d.σ
decay(d::EGPpowermix) = d.ξ    # noms à vérifier...
shape_tail(d::EGPpowermix) = d.κ₁
shape_central(d::EGPpowermix) = d.κ₂
prob(d::EGPpowermix) = d.p

params(d::EGPpowermix) = (d.σ, d.ξ, d.κ₁, d.κ₂, d.p)
partype(::EGPpowermix{T}) where {T} = T


#### Evaluation

function logpdf(d::EGPpowermix{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, κ₁, κ₂, p) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    g(v::Real, p::Real, κ₁::Real, κ₂::Real) = p*κ₁*v^(κ₁-1) + (1-p)*κ₂*v^(κ₂-1)
    lg(v::Real, p::Real, κ₁::Real, κ₂::Real) = log(g(v,p,κ₁,κ₂))

    lf = -log(σ) + Distributions.logpdf(pd, x/σ) + lg(Distributions.cdf(pd, x/σ),p,κ₁,κ₂)

    return lf
end

pdf(d::EGPpowermix, x::Real) = exp(logpdf(d, x))

function logcdf(d::EGPpowermix{T}, x::Real) where T<:Real
    μ = 0
    (σ, ξ, κ₁, κ₂, p) = params(d)

    pd = GeneralizedPareto(μ, 1, ξ)

    G(v::Real, p::Real, κ₁::Real, κ₂::Real) = p*v^κ₁ + (1-p)*v^κ₂

    lf = log(G(Distributions.cdf(pd, x/σ), p, κ₁, κ₂))

    return lf
end

cdf(d::EGPpowermix, x::Real) = exp(logcdf(d, x))

function quantile(d::EGPpowermix{T}, p::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    (σ, ξ, κ₁, κ₂, m) = params(d)

    G(v::Real, m::Real, κ₁::Real, κ₂::Real) = m*v^κ₁ + (1-m)*v^κ₂

    function invG(p, m, κ₁, κ₂)
        a = 0
        b = 1
        δ = 0.00001
        crit = 1
        X = NaN
        while crit > 2*δ
            X = (a + b)/2
            F = G(X, m, κ₁, κ₂)
            if F <= p
                a = X
            else
                b = X
            end
            crit = b - a
        end
        return X
    end

    x = (σ/ξ)*(((1 - invG(p,m,κ₁, κ₂))^(-ξ)) - 1)

    return x
end


#### Sampling

function rand(rng::AbstractRNG, d::EGPpowermix)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    G(v::Real, m::Real, κ₁::Real, κ₂::Real) = m*v^κ₁ + (1-m)*v^κ₂

    function invG(p, m, κ₁, κ₂)
        a = 0
        b = 1
        δ = 0.00001
        crit = 1
        X = NaN
        while crit > 2*δ
            X = (a + b)/2
            F = G(X, m, κ₁, κ₂)
            if F <= p
                a = X
            else
                b = X
            end
            crit = b - a
        end
        return X
    end

    return (d.σ/d.ξ)*(((1 - invG(u,d.p,d.κ₁,d.κ₂))^(-d.ξ)) - 1)
end

sampler(d::EGPpowermix) = d
