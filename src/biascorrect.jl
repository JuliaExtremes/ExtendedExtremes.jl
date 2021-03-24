"""
    T_transform(obs_ref, model_ref, model_fut, x) (stationnary)

"""
function T_transform(obs_ref::ContinuousUnivariateDistribution,  # dist des obs période de référence
                     model_ref::ContinuousUnivariateDistribution,  # dist du modèle période de référence
                     model_fut::ContinuousUnivariateDistribution,  # dist du modèle période future
                     x::Real)

    Fmodel_fut = cdf(model_fut, x)

    Fm1model_ref = quantile(model_ref, Fmodel_fut)

    Fobs = cdf(obs_ref,  Fm1model_ref)

    return Fobs
end

"""
    T_transform(obs_ref, model_ref, model_fut, covar_ref, covar_fut, x)  (non-stationnary)

"""
function T_transform(obs_ref::nonstatEGPpower,  # dist des obs période de référence
                     model_ref::nonstatEGPpower,  # dist du modèle période de référence
                     model_fut::nonstatEGPpower,  # dist du modèle période future
                     covar_ref::Real,  # covariable pour la période de référence (ex : [CO₂] de l'année 1980)
                     covar_fut::Real, # covariable pour la période future (ex : [CO₂] de l'année 2050)
                     x::Real)

    Fmodel_fut = cdf(model_fut, x, covar_fut)

    Fm1model_ref = quantile(model_ref, Fmodel_fut, covar_ref)

    Fobs = cdf(obs_ref,  Fm1model_ref, covar_ref)

    return Fobs
end

"""
    invT_transform(obs_ref, model_ref, model_fut, p) (stationnary)

"""
function invT_transform(obs_ref::ContinuousUnivariateDistribution,  # dist des obs période de référence
                        model_ref::ContinuousUnivariateDistribution,  # dist du modèle période de référence
                        model_fut::ContinuousUnivariateDistribution,  # dist du modèle période future
                        p::Real)

    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    # Méthode de la bisection

    a = 0
    b = 1000  # à vérifier
    δ = 0.00001
    crit = 1
    X = NaN
    while crit > 2*δ
        X = (a + b)/2
        F = T_transform(obs_ref, model_ref, model_fut, X)
        if F <= p
            a = X
        else
            b = X
        end
        crit = b - a
    end
    return X
end

"""
    invT_transform(obs_ref, model_ref, model_fut, covar_ref, covar_fut, p)  (non-stationnary)

"""
function invT_transform(obs_ref::nonstatEGPpower,  # dist des obs période de référence
              model_ref::nonstatEGPpower,  # dist du modèle période de référence
              model_fut::nonstatEGPpower,  # dist du modèle période future
              covar_ref::Real,  # covariable pour la période de référence (ex : [CO₂] de l'année 1980)
              covar_fut::Real, # covariable pour la période future (ex : [CO₂] de l'année 2050)
              p::Real)

    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."


    # Méthode de la bisection

    a = 0
    b = 1000  # à vérifier
    δ = 0.00001
    crit = 1
    X = NaN
    while crit > 2*δ
        X = (a + b)/2
        F = T_transform(obs_ref, model_ref, model_fut, covar_ref, covar_fut, X)
        if F <= p
            a = X
        else
            b = X
        end
        crit = b - a
    end
    return X
end


# Développement de l'algo de correction :
# On a besoin de obs_ref, model_ref et model_fut, on va commencer par le cas stationnaire
# Étapes :
# 1. Ajuster le modèle sur les observations pour obtenir obs_ref
# 2. Ajuster le modèle sur les données de calibration du modèle pour obtenir model_ref
# 3. Ajuster le modèle sur les données futures du modèle pour obtenir model_fut
# 4. Calculer la cdf Fmodel_fut
# 5. Calculer la cdf-transform inverse de Fmodel_fut avec obs_ref, model_ref et model_fut
# 6. Comparer les valeurs obtenues avec la valeurs simulées
