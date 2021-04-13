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
