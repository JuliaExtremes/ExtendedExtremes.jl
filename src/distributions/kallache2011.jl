"""
    CDFT(Yc, Xc, Xp)

"""
struct CDFT
    Yc::ContinuousUnivariateDistribution  # Distribution des données observées pour la période de calibration
    Xc::ContinuousUnivariateDistribution  # Distribution des données du modèle pour la période de calibration
    Xp::ContinuousUnivariateDistribution  # Distribution des données du modèle pour la période projetée
end


"""
    T_transform(m::CDFT, x::Real)

"""
function T_transform(m::CDFT, x::Real)

    F_Xp = cdf(m.Xp, x)  # Évaluation en x de la CDF des données du modèle pour la période projeté

    Fm1_Xc = quantile(m.Xc, F_Xp)  # Quantile correspondant de la distribution des données du modèle pour la période de calibration

    F_Yp = cdf(m.Yc,  Fm1_Xc)  # Évaluation de la CDF des données observées pour la période de calibration au point correspondant

    return F_Yp  # Valeur de la CDF évaluée au point x
end

#### Evaluation

cdf(m::CDFT, x::Real) = T_transform(m, x)

function pdf(m::CDFT, x::Real)
    # très approximatif... à revoir
    # Inf en zéro (?)

    δ = 0.00001
    x₁ = x - δ
    x₂ = x + δ

    return (cdf(m, x₂) - cdf(m, x₁)) / (x₂ - x₁)
end

function quantile(m::CDFT,
                  p::Real;
                  interval::Array{Float64,1}=[0.0, 500.0])

    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    res = optimize(θ -> (p - cdf(m, θ))^2, interval[1], interval[2])

    return Optim.minimizer(res)[1]
end

#### Sampling

function rand(rng::AbstractRNG, m::CDFT)
    # Generate a Float64 random number uniformly in (0,1].
    u = 1 - rand(rng)

    return quantile(m, u)
end

sampler(m::CDFT) = m    # à vérifier
