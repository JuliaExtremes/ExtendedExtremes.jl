# TO-DO :
# -tenir compte du nombre de parametres pour les initialvalues dans la fct fit
# -vérifier la convergence pour l'optimisation dans la fct fit


"""
    nonstatEGPpower(σ, ξ₀, ξ₁, κ₀, κ₁)

Extension non stationnaire du modèle 1 de Naveau. À détailler...
"""
struct nonstatEGPpower{T<:Real} <: ContinuousUnivariateDistribution
    σ::T    # scale parameter
    ξ₀::T    # rate of upper tail decay
    ξ₁::T    # rate of upper tail decay
    κ₀::T    # shape of the lower tail
    κ₁::T    # shape of the lower tail

    function nonstatEGPpower{T}(σ::T, ξ₀::T, ξ₁::T, κ₀::T, κ₁::T) where {T <: Real}
        new{T}(σ, ξ₀, ξ₁, κ₀, κ₁)
    end
end

function nonstatEGPpower(σ::T, ξ₀::T, ξ₁::T, κ₀::T, κ₁::T; check_args=true) where {T <: Real}
    #check_args && @check_args(EGPpower, σ > zero(σ) && κ > zero(κ))
    return nonstatEGPpower{T}(σ, ξ₀, ξ₁, κ₀, κ₁)
end

nonstatEGPpower(σ::Real, ξ₀::Real, ξ₁::Real, κ₀::Real, κ₁::Real) = nonstatEGPpower(promote(σ, ξ₀, ξ₁, κ₀, κ₁)...)
nonstatEGPpower(σ::Integer, ξ₀::Integer, ξ₁::Integer, κ₀::Integer, κ₁::Integer) = nonstatEGPpower(float(σ), float(ξ₀), float(ξ₁), float(κ₀), float(κ₁))

params(d::nonstatEGPpower) = (d.σ, d.ξ₀, d.ξ₁, d.κ₀, d.κ₁)
#### Evaluation

function logpdf(d::nonstatEGPpower{T}, x::Real, covariate::Real) where T<:Real
    μ = 0
    σ = d.σ
    ξ = d.ξ₀ + d.ξ₁*covariate
    κ = d.κ₀ + d.κ₁*covariate

    pd = GeneralizedPareto(μ, 1, ξ)

    lg(v::Real) = log(κ) + (κ-1)*log(v)

    p = -log(σ) + Distributions.logpdf(pd, x/σ) + lg(Distributions.cdf(pd, x/σ))

    return p
end

pdf(d::nonstatEGPpower, x::Real, covariate::Real) = exp(logpdf(d, x, covariate))

function logcdf(d::nonstatEGPpower{T}, x::Real, covariate::Real) where T<:Real
    μ = 0
    σ = d.σ
    ξ = d.ξ₀ + d.ξ₁*covariate
    κ = d.κ₀ + d.κ₁*covariate

    pd = GeneralizedPareto(μ, 1, ξ)

    #G(v::Real, κ::Real) = v^κ
    #lG(v::Real, κ::Real) = κ*log(v)

    p = κ*Distributions.logcdf(pd, x/σ)

    return p
end

cdf(d::nonstatEGPpower, x::Real, covariate::Real) = exp(logcdf(d, x, covariate))

function quantile(d::nonstatEGPpower{T}, p::Real, covariate::Real) where T<:Real
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    σ = d.σ
    ξ = d.ξ₀ + d.ξ₁*covariate
    κ = d.κ₀ + d.κ₁*covariate

    invG(p::Real, κ::Real) = p^(1/κ)

    x = (σ/ξ)*(((1 - invG(p,κ))^(-ξ)) - 1)

    return x
end

#### Parameter estimation

function EGPpowerfit(data::Array{<:Real,1},
    covariate::Array{<:Real,1};  # Vecteur de la variable explicative
    initialvalues::Vector{<:Real}=Float64[],
    lowertailcov::Bool=true,  # variable explicative sur κ ?
    uppertailcov::Bool=true,  # variable explicative sur ξ ?
    censoring::Real=0)

    # TO-DO :
    # -tenir compte du nombre de parametres pour les initialvalues
    # -voir pour ajouter la censure

    if isempty(initialvalues)
        σ₀, ξ₀₀, ξ₁₀, κ₀₀, κ₁₀ = [1, 0.15, 0, 1, 1]
    else
        σ₀, ξ₀₀, ξ₁₀, κ₀₀, κ₁₀ = initialvalues
    end

    if !lowertailcov & !uppertailcov
        return EGPpowerfit(data, initialvalues = [σ₀, ξ₀₀, κ₀₀], censoring = censoring)
    else
        r_data = data[data .>= censoring]
        l_data = fill(censoring, count(data .< censoring))

        r_cov = covariate[data .>= censoring]
        l_cov = covariate[data .< censoring]

        # Modèle avec ξ = ξ₀ + ξ₁*covariate
        function loglike_xi(σ, ξ₀, ξ₁, κ₀)
            if σ <= 0 || κ₀ <= 0
                return -Inf
            else
                pd = nonstatEGPpower(σ, ξ₀, ξ₁, κ₀, 0)
                #ll = sum(logpdf.(pd, data, covariate))
                ll = sum(logcdf.(pd, l_data, l_cov)) + sum(logpdf.(pd, r_data, r_cov))
                return ll
            end
        end
        # Modèle avec κ = κ₀ + κ₁*covariate
        function loglike_kappa(σ, ξ₀, κ₀, κ₁)
            if σ <= 0 || κ₀*κ₁ <= 0
                return -Inf
            else
                pd = nonstatEGPpower(σ, ξ₀, 0, κ₀, κ₁)
                #ll = sum(logpdf.(pd, data, covariate))
                ll = sum(logcdf.(pd, l_data, l_cov)) + sum(logpdf.(pd, r_data, r_cov))
                return ll
            end
        end
        # Modèle avec ξ = ξ₀ + ξ₁*covariate et κ = κ₀ + κ₁*covariate
        function loglike_all(σ, ξ₀, ξ₁, κ₀, κ₁)
            if σ <= 0 || κ₀*κ₁ <= 0
                return -Inf
            else
                pd = nonstatEGPpower(σ, ξ₀, ξ₁, κ₀, κ₁)
                #ll = sum(logpdf.(pd, data, covariate))
                ll = sum(logcdf.(pd, l_data, l_cov)) + sum(logpdf.(pd, r_data, r_cov))
                return ll
            end
        end

        # Fonctions objectives
        fobj_1(θ) = -loglike_xi(θ...)
        fobj_2(θ) = -loglike_kappa(θ...)
        fobj_3(θ) = -loglike_all(θ...)

        if !lowertailcov & uppertailcov  # Modèle avec ξ = ξ₀ + ξ₁*covariate
            #println("Model 1")
            res = optimize(fobj_1, [σ₀, ξ₀₀, ξ₁₀, κ₀₀])
            fittedmodel = nonstatEGPpower(Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3], Optim.minimizer(res)[4], 0)
        elseif lowertailcov & !uppertailcov  # Modèle avec κ = κ₀ + κ₁*covariate
            #println("Model 2")
            res = optimize(fobj_2, [σ₀, ξ₀₀, κ₀₀, κ₁₀])
            fittedmodel = nonstatEGPpower(Optim.minimizer(res)[1], Optim.minimizer(res)[2], 0, Optim.minimizer(res)[3], Optim.minimizer(res)[4])
        else  # Modèle avec ξ = ξ₀ + ξ₁*covariate et κ = κ₀ + κ₁*covariate
            #println("Model 3")
            res = optimize(fobj_3, [σ₀, ξ₀₀, ξ₁₀, κ₀₀, κ₁₀])
            fittedmodel = nonstatEGPpower(Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3], Optim.minimizer(res)[4], Optim.minimizer(res)[5])
        end

        return fittedmodel, getparametervalues.(fittedmodel, covariate)

        # TO-DO : - vérifier la convergence, ex :
            #if Optim.converged(res)
            #    σ̂, ξ̂, κ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3]]
            #else
            #    @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
            #    σ̂, ξ̂, κ̂  = [initialvalues[1], initialvalues[2], initialvalues[3]]
            #end
    end
end


"""
     getparametervalues(fm, covariate)

Pour obtenir les valeurs des paramètres à partir de la valeur de la variable explicative.
"""
function getparametervalues(fm::nonstatEGPpower{T}, covariate::Real) where T<:Real
    σ = fm.σ
    ξ = fm.ξ₀ + (fm.ξ₁ * covariate)
    κ = fm.κ₀ + (fm.κ₁ * covariate)
    return EGPpower(σ, ξ, κ)
end


"""
    BIC(fm, data, covariate)

"""
function BIC(fm::nonstatEGPpower{T}, data::Array{T,1}, covariate::Array{T,1}; censoring::Real=0) where T<:Real
    n = size(data,1)
    k = size(params(fm),1)
    fm.ξ₁ == 0 ? k = k-1 : nothing
    fm.κ₁ == 0 ? k = k-1 : nothing

    r_data = data[data .>= censoring]
    l_data = fill(censoring, count(data .< censoring))

    r_cov = covariate[data .>= censoring]
    l_cov = covariate[data .< censoring]

    ll = sum(logcdf.(fm, l_data, l_cov)) + sum(logpdf.(fm, r_data, r_cov))

    return ll - k/2*log(n)
end
