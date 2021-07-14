"""
    T_transform(Yc, Xc, Xp, x) (stationnary)

"""
function T_transform(Yc::ContinuousUnivariateDistribution,  # Distribution des données observées pour la période de calibration
                     Xc::ContinuousUnivariateDistribution,  # Distribution des données du modèle pour la période de calibration
                     Xp::ContinuousUnivariateDistribution,  # Distribution des données du modèle pour la période projetée
                     x::Real)

    F_Xp = cdf(Xp, x)  # Évaluation en x de la CDF des données du modèle pour la période projeté

    Fm1_Xc = quantile(Xc, F_Xp)  # Quantile correspondant de la distribution des données du modèle pour la période de calibration

    F_Yp = cdf(Yc,  Fm1_Xc)  # Évaluation de la CDF des données observées pour la période de calibration au point correspondant

    return F_Yp  # Valeur de la CDF évaluée au point x
end


"""
    T_transform(Yx, Xc, Xc_cov, Xp, Xp_cov, x)  (non-stationnary)

"""
function T_transform(Yc::EGPpower,  # Distribution des données observées pour la période de calibration
                     Xc::nonstatEGPpower,  # Distribution des données du modèle pour la période de calibration
                     Xc_cov::Real,  # covariable pour la période de calibration (ex : [CO₂] de l'année 1980)
                     Xp::nonstatEGPpower,  # Distribution des données du modèle pour la période projetée
                     Xp_cov::Real,  # covariable pour la période projetée (ex : [CO₂] de l'année 2050)
                     x::Real)

    F_Xp = cdf(Xp, x, Xp_cov)  # Évaluation en x de la CDF des données du modèle pour la période projeté

    Fm1_Xc = quantile(Xc, F_Xp, Xc_cov)  # Quantile correspondant de la distribution des données du modèle pour la période de calibration

    F_Yp = cdf(Yc,  Fm1_Xc)  # Évaluation de la CDF des données observées pour la période de calibration au point correspondant

    return F_Yp  # Valeur de la CDF évaluée au point x
end

## Bias correction

"""
    bias_correct(obs, ref, covariate_ref; fut, covariate_fut, censoring_obs, censoring_sim, wet_thresh)

"""
function bias_correct(obs::Array{Float64, 1},
                ref::Array{Float64, 1},
                covariate_ref::Array{Float64, 1};
                fut::Array{Float64, 1}=Float64[],
                covariate_fut::Array{Float64, 1}=Float64[],
                censoring_obs::Real=5.0,
                censoring_sim::Real=5.0,
                wet_thresh::Float64=1.0)

    # Transformation des valeurs en fonction du seuil
    obs_t = obs[obs .> wet_thresh] .- wet_thresh

    ref_t = ref[ref .> wet_thresh] .- wet_thresh
    covariate_ref_t = covariate_ref[ref .> wet_thresh]

    fut_t = fut[fut .> wet_thresh] .- wet_thresh
    covariate_fut_t = covariate_fut[fut .> wet_thresh]

    # Ajustement du modèle sur les observations
    fm_obs = EGPpowerfit(obs_t, censoring=censoring_obs)

    # Fusion des prédictions et des projections pour l'ajustement
    sim = vcat(ref_t, fut_t)
    covariate_sim = vcat(covariate_ref_t, covariate_fut_t)

    # Ajustement du modèle sur les simulations
    fm_sim, fms = EGPpowerfit(sim, censoring=censoring_sim, covariate=covariate_sim, scalecov=false, uppertailcov=true, lowertailcov=false)

    # Correction des prédictions
    output_ref = similar(ref)
    output_ref[ref .<= wet_thresh] .= 0.0
    output_ref[ref .> wet_thresh] = quantile.(fm_obs, cdf.(EGPpower.(fm_sim, covariate_ref_t), ref_t)) .+ wet_thresh

    if !isempty(fut)
        # Récupération de la CDF des obs pour le futur (Kallache)
        Yp = CDFT.(Ref(fm_obs), Ref(EGPpower(fm_sim, last(covariate_ref))), EGPpower.(fm_sim, covariate_fut_t))
        p = cdf.(getfield.(Yp, :Xp), fut_t)

        # Correction des projections
        output_fut = similar(fut)
        output_fut[fut .<= wet_thresh] .= 0.0
        output_fut[fut .> wet_thresh] = quantile.(Yp, p) .+ wet_thresh
        append!(output_ref, output_fut)

    end
    return output_ref
end

function bias_correct(obs::Array{Float64, 1},
                ref::Array{Float64, 1},
                covariate_ref::Array{Float64, 1},
                fut::Array{Float64, 1},
                covariate_fut::Array{Float64, 1};
                censoring_obs::Real=5.0,
                censoring_sim::Real=5.0,
                wet_thresh::Float64=1.0)

    return bias_correct(obs, ref, covariate_ref, fut=fut, covariate_fut=covariate_fut, censoring_obs=censoring_obs, censoring_sim=censoring_sim, wet_thresh=wet_thresh)
end

function bias_correct(obs::TimeArray{Float64,1,T,Array{Float64,1}},
                ref::TimeArray{Float64,1,T,Array{Float64,1}},
                covariate_ref::TimeArray{Float64,1,T,Array{Float64,1}},
                fut::TimeArray{Float64,1,T,Array{Float64,1}},
                covariate_fut::TimeArray{Float64,1,T,Array{Float64,1}};
                censoring_obs::Real=5.0,
                censoring_sim::Real=5.0,
                wet_thresh::Float64=1.0) where {T <: TimeType}

    x = bias_correct(values(obs), values(ref), values(covariate_ref), fut=values(fut), covariate_fut=values(covariate_fut), censoring_obs=censoring_obs,
      censoring_sim=censoring_sim, wet_thresh=wet_thresh)
    t = vcat(timestamp(ref), timestamp(fut))

    return TimeArray(t, x)
end

function bias_correct(obs::TimeArray{Float64,1,T,Array{Float64,1}},
                ref::TimeArray{Float64,1,T,Array{Float64,1}},
                covariate_ref::TimeArray{Float64,1,T,Array{Float64,1}};
                censoring_obs::Real=5.0,
                censoring_sim::Real=5.0,
                wet_thresh::Float64=1.0) where {T <: TimeType}

    x = bias_correct(values(obs), values(ref), values(covariate_ref), censoring_obs=censoring_obs, censoring_sim=censoring_sim, wet_thresh=wet_thresh)
    t = timestamp(ref)

    return TimeArray(t, x)
end
