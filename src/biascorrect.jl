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





# """
#     bias_correct(obs_ts, obs_timevec, sim_ts, sim_timevec, sim_covariate)
#
# """
# function bias_correct(obs_ts::Array{Float64,1},  # Série temporelle des observations
#                       obs_timevec::Array{Int64,1},  # Vecteur des années associé
#                       sim_ts::Array{Float64,1},  # Série temporelle des valeurs simulées
#                       sim_timevec::Array{Int64,1},  # Vecteur des années associé
#                       sim_covariate::Array{Float64,1};  # Vecteur des covariables ([CO₂])
#                       censoring_obs::Real=2,
#                       censoring_sim::Real=1)
#
#     # Quelques vérifications de dimensions pour les série temporelle/vecteurs d'années
#     @assert length(obs_ts) == length(obs_timevec)
#     @assert length(sim_ts) == length(sim_timevec)
#     @assert length(sim_ts) == length(sim_covariate)
#
#     # On récupère les années
#     years_obs = unique(obs_timevec)
#     years_sim = unique(sim_timevec)
#
#
#     # Ajustement des modèles sur les données :
#     fm_obs = EGPpowerfit(obs_ts, censoring=censoring_obs)  # Ajustement du modèle sur obs (la censure sera un hyperparamètre éventuellement)
#     fm_sim, fms_sim = EGPpowerfit(sim_ts, sim_covariate, censoring=censoring_sim);  # Ajustement du modèle ns sur simu (la censure sera un hyperparamètre éventuellement)
#
#     # Boucle de correction
#     corrected_ts = Array{Float64}(undef, length(sim_ts))  # Génération d'un vecteur vide de la dimension des valeurs simulées
#     k=0
#     for y in years_sim  # Pour chaque année simulée,
#
#         # On regarde si cette année existe dans les observations
#
#         if y in years_obs  # Si oui,
#
#             # On récupère les valeurs simulées pour cette année
#             idx = sim_timevec .== y
#             ts_y = sim_ts[idx]
#             co2_y = unique(sim_covariate[idx])[1]  # ainsi que la [CO₂]
#
#             cs = Float64[]  # Vecteur vide pour stocker les valeurs corrigées
#
#             for i in ts_y  # Pour chaque valeur simulée,
#
#                 cv = cdf(fm_sim, i, co2_y)  # On évalue la probabilité cumulative du modèle en cette valeur
#
#                 # On applique Kallache (on cherche la valeur qui fait que les probabilités cumulatives des deux modèles sont égales)
#                 f(x::Float64) = (cv - T_transform(fm_obs, fm_sim, co2_y, fm_sim, co2_y, x))^2
#
#                 # Cela revient à trouver la valeur qui minimise l'écart
#                 res = optimize(f, 0.0, 250.0)
#
#                 # On ajoute cette valeur au vecteur des valeurs corrigées
#                 push!(cs, Optim.minimizer(res)[1])
#             end
#
#             # On vient ajouter toutes les valeurs pour cette année au grand vecteurs des valeurs simulées
#             n = length(cs)
#             for j = 1:n
#                 corrected_ts[k+j] = cs[j]
#             end
#             k += n
#
#         elseif y > years_obs[end]  # Si l'année n'existe dans les observations, on est dans le "futur"
#
#             # On récupère les valeurs simulées pour cette année
#             idx = sim_timevec .== y
#             ts_y = sim_ts[idx]
#             co2_y = unique(sim_covariate[idx])[1]
#
#             # Et les valeurs pour la dernière année observée
#             idx = sim_timevec .== years_obs[end]
#             ts_e = sim_ts[idx]
#             co2_e = unique(sim_covariate[idx])[1]
#
#             cs = Float64[]  # Vecteur vide pour stocker les valeurs corrigées
#
#             for i in ts_y  # Pour chaque valeur simulée,
#
#                 cv = cdf(fm_sim, i, co2_y)  # On évalue la probabilité cumulative du modèle en cette valeur
#
#                 # On applique Kallache (on cherche la valeur qui fait que les probabilités cumulatives des deux modèles sont égales)
#                 f(x::Float64) = (cv - T_transform(fm_obs, fm_sim, co2_e, fm_sim, co2_y, x))^2
#
#                 # Cela revient à trouver la valeur qui minimise l'écart
#                 res = optimize(f, 0.0, 250.0)
#
#                 # On ajoute cette valeur au vecteur des valeurs corrigées
#                 push!(cs, Optim.minimizer(res)[1])
#             end
#
#             # On vient ajouter toutes les valeurs pour cette année au grand vecteurs des valeurs simulées
#             n = length(cs)
#             for j = 1:n
#                 corrected_ts[k+j] = cs[j]
#             end
#             k += n
#         end
#     end
#     return corrected_ts  # On retourne le vecteur des valeurs corrigées
# end


"""
    bias_correct(obs, sim, sim_covariate)

"""
function bias_correct(obs::TimeArray{Float64,1,T,Array{Float64,1}},  # Série temporelle des observations
                      sim::TimeArray{Float64,1,T,Array{Float64,1}}, # Série temporelle des valeurs simulées
                      sim_covariate::TimeArray{Float64,1,T,Array{Float64,1}};  # Vecteur des covariables ([CO₂])
                      censoring_obs::Real=1,
                      censoring_sim::Real=1) where {T <: TimeType}

    # On récupère les années
    years_obs = unique(year.(timestamp(obs)))
    years_sim = unique(year.(timestamp(sim)))

    # Ajustement des modèles sur les données :
    fm_obs = EGPpowerfit(values(obs), censoring=censoring_obs)  # Ajustement du modèle sur obs (la censure sera un hyperparamètre éventuellement)
    fm_sim, fms_sim = EGPpowerfit(values(sim), values(sim_covariate), censoring=censoring_sim);  # Ajustement du modèle ns sur simu (la censure sera un hyperparamètre éventuellement)
    fm_sim = collapse(TimeArray(timestamp(sim), fms_sim), year, first)

    # KALLACHE
    Yc = fill(fm_obs, length(years_sim))
    Xp = values(fm_sim)
    Xc = copy(Xp)
    Xc[years_sim .> last(years_obs)] .= values(when(fm_sim, year, last(years_obs)))[1]
    Yp = CDFT.(Yc, Xc, Xp)

    # CORRECTION
    corrected_values = Float64[]
    sizehint!(corrected_values, length(sim))
    for i in eachindex(years_sim)
        x = when(sim, year, years_sim[i]);
        p = cdf.(Yp[i].Xp, values(x));
        ŷ = quantile.(Yp[i], p)
        append!(corrected_values, ŷ)
    end

    return TimeArray(timestamp(sim), corrected_values)
end
