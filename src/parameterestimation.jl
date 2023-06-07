
function fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y::Vector{<:Real}, initialvalues::Vector{<:Real}; leftcensoring::Real)
   
    # Values above the censoring threshold
    y⁺ = filter( v -> v > leftcensoring, y)

    # Number of values below the censoring threshold
    n⁻ = count( y .< leftcensoring)
    
    V = EGPtype(pd)
    
    function loglike(κ::Real, σ::Real, ξ::Real)
        if (κ > 0) & (σ > 0)
            pd = ExtendedGeneralizedPareto(V(κ), GeneralizedPareto(σ, ξ))
#             return sum(logpdf.(pd, y)) - (κ - 1.)^2/.1
            if n⁻ == 0
                return sum(logpdf.(pd, y⁺)) - (κ - 1.)^2/.1
            else
                return sum(logpdf.(pd, y⁺)) - (κ - 1.)^2/.1 + n⁻ * logcdf(pd, leftcensoring)
            end
        else
            return -Inf
        end
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        κ̂, σ̂, ξ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3]]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        κ̂, σ̂, ξ̂   = [initialvalues[1], initialvalues[2], initialvalues[3]]
    end
    
    return ExtendedGeneralizedPareto(V(κ̂), GeneralizedPareto(σ̂, ξ̂))
    
end

function fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y::Vector{<:Real}; leftcensoring::Real=0.)
   
    initialvalues = [1.0, 1.0, 0.0]
    
    return fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y, initialvalues, leftcensoring = leftcensoring)
    
end





# function fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y::Vector{<:Real}, initialvalues::Vector{<:Real})
   
#     #κ₀, σ₀, ξ₀, = initialvalues
    
#     V = EGPtype(pd)
    
#     function loglike(κ::Real, σ::Real, ξ::Real)
#         if (κ > 0) & (σ > 0)
#             pd = ExtendedGeneralizedPareto(V(κ), GeneralizedPareto(σ, ξ))
#             return sum(logpdf.(pd, y))
#         else
#             return -Inf
#         end
#     end

#     fobj(θ) = -loglike(θ...)

#     res = optimize(fobj, initialvalues)

#     if Optim.converged(res)
#         κ̂, σ̂, ξ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3]]
#     else
#         @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#         κ̂, σ̂, ξ̂   = [initialvalues[1], initialvalues[2], initialvalues[3]]
#     end
    
#     return ExtendedGeneralizedPareto(V(κ̂), GeneralizedPareto(σ̂, ξ̂))
    
# end

# function fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y::Vector{<:Real})
   
#     initialvalues = [1.0, 1.0, 0.0]
    
#     return fit_mle(pd::Type{<:ExtendedGeneralizedPareto}, y, initialvalues)
    
# end

# """
#     EGPpowerfit()

# Estimate the first extended GP model of Naveau et al. (2016) parameters by maximum likelihood.
# """
# function EGPpowerfit(data::Array{<:Real,1};
#     initialvalues::Vector{<:Real}=Float64[],
#     censoring::Real=0,
#     covariate::Array{<:Real,1}=Float64[],  # Vecteur de la variable explicative
#     scalecov::Bool=true, # variable explicative sur σ ?
#     lowertailcov::Bool=false,  # variable explicative sur κ ?
#     uppertailcov::Bool=false)  # variable explicative sur ξ ?)

#     if !isempty(covariate)
#         return EGPnonstatpowerfit(data, covariate, initialvalues=initialvalues, scalecov=scalecov, lowertailcov=lowertailcov, uppertailcov=uppertailcov, censoring=censoring)
#     else
#         if isempty(initialvalues)
#             σ₀, ξ₀, κ₀ = [1, 0.15, 1]
#             initialvalues = [σ₀, ξ₀, κ₀]
#         end

#         function loglike(σ::Real, ξ::Real, κ::Real)
#             if κ <= 0 || σ <= 0
#                 return -Inf
#             else
#                 pd = EGPpower(σ, ξ, κ)

#                 r_data = data[data .>= censoring]
#                 l_data = fill(censoring, count(data .< censoring))

#                 ll = sum(logcdf.(pd, l_data)) + sum(logpdf.(pd, r_data))

#                 return ll
#             end
#         end

#         fobj(θ) = -loglike(θ...)

#         res = optimize(fobj, initialvalues)

#         if Optim.converged(res)
#             σ̂, ξ̂, κ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3]]
#         else
#             @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#             σ̂, ξ̂, κ̂  = [initialvalues[1], initialvalues[2], initialvalues[3]]
#         end

#         fittedmodel = EGPpower(σ̂, ξ̂, κ̂)

#         return fittedmodel
#     end
# end

# function EGPpowerfit(data::TimeArray;
#     initialvalues::Vector{<:Real}=Float64[],
#     censoring::Real=0,
#     covariate::TimeArray=TimeArray(DateTime[], Float64[]),  # Vecteur de la variable explicative
#     scalecov::Bool=true,
#     lowertailcov::Bool=false,  # variable explicative sur κ ?
#     uppertailcov::Bool=false)  # variable explicative sur ξ ?)

#     if isempty(covariate)
#         return EGPpowerfit(values(data), initialvalues=initialvalues, scalecov=scalecov, lowertailcov=lowertailcov, uppertailcov=uppertailcov, censoring=censoring)
#     else
#         return EGPpowerfit(values(data), covariate=values(covariate), initialvalues=initialvalues,
#         scalecov=scalecov, lowertailcov=lowertailcov, uppertailcov=uppertailcov, censoring=censoring)
#     end
# end

# """
#     EGPpowermixfit()

# Estimate the second extended GP model of Naveau et al. (2016) parameters by maximum likelihood.
# """
# function EGPpowermixfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)

#     #inner_optimizer = NelderMead()
#     #lower = [0, -Inf, 0, 0, 0]
#     #upper = [Inf, Inf, Inf, Inf, 1]

#     if isempty(initialvalues)
#         σ₀, ξ₀, κ₁₀, κ₂₀, p₀ = [1, .15, 1, 2, .4]  # à remplacer par une fonction (getinitialvalue)
#         initialvalues = [σ₀, ξ₀, κ₁₀, κ₂₀, p₀]
#     end

#     function loglike(σ::Real, ξ::Real, κ₁::Real, κ₂::Real, p::Real)
#         if p < 0 || p > 1 || κ₁ <= 0 || κ₂ <= 0 || σ <= 0 || κ₂ < κ₁
#             return -Inf
#         else
#             pd = EGPpowermix(σ, ξ, κ₁, κ₂, p)

#             r_data = data[data .>= censoring]
#             l_data = fill(censoring, count(data .< censoring))

#             ll = sum(logcdf.(pd, l_data)) + sum(logpdf.(pd, r_data))

#             return ll
#         end
#     end

#     fobj(θ) = -loglike(θ...)

#     res = optimize(fobj, initialvalues)

#     if Optim.converged(res)
#         σ̂, ξ̂, κ̂₁, κ̂₂, p̂ = [Optim.minimizer(res)...]
#     else
#         @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#         σ̂, ξ̂, κ̂₁, κ̂₂, p̂  = [initialvalues...]
#     end

#     fittedmodel = EGPpowermix(σ̂, ξ̂, κ̂₁, κ̂₂, p̂)

#     return fittedmodel
# end

# function EGPpowermixfit(data::TimeArray; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)
#     return EGPpowermixfit(values(data), initialvalues=initialvalues, censoring=censoring)
# end


# """
#     EGPbetafit()

# Estimate the third extended GP model of Naveau et al. (2016) parameters by maximum likelihood.
# """
# function EGPbetafit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)

#     if isempty(initialvalues)
#         σ₀, ξ₀, δ₀ = [5, 0.15, 2]
#         initialvalues = [σ₀, ξ₀, δ₀]
#     end

#     function loglike(σ::Real, ξ::Real, δ::Real)
#         if δ <= 0 || δ > 100 || σ <= 0 || ξ > 0.99
#             return -Inf
#         else
#             pd = EGPbeta(σ, ξ, δ)

#             r_data = data[data .>= censoring]
#             l_data = fill(censoring, count(data .< censoring))

#             ll = sum(logcdf.(pd, l_data)) + sum(logpdf.(pd, r_data))

#             return ll
#         end
#     end

#     fobj(θ) = -loglike(θ...)

#     res = optimize(fobj, initialvalues)

#     if Optim.converged(res)
#         σ̂, ξ̂, δ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3]]
#     else
#         @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#         σ̂, ξ̂, δ̂  = [initialvalues[1], initialvalues[2], initialvalues[3]]
#     end

#     fittedmodel = EGPbeta(σ̂, ξ̂, δ̂)

#     return fittedmodel
# end

# function EGPbetafit(data::TimeArray; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)
#     return EGPbetafit(values(data), initialvalues=initialvalues, censoring=censoring)
# end


# """
#     EGPbetapowerfit()

# Estimate the fourth extended GP model of Naveau et al. (2016) parameters by maximum likelihood.
# """
# function EGPbetapowerfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)

#     if isempty(initialvalues)
#         σ₀, ξ₀, δ₀, κ₀ = [1, 0.15, 1, 1]
#         initialvalues = [σ₀, ξ₀, δ₀, κ₀]
#     end

#     function loglike(σ::Real, ξ::Real, δ::Real, κ::Real)
#         if κ <= 0 || δ <= 0 || δ > 100 || σ <= 0
#             return -Inf
#         else
#             pd = EGPbetapower(σ, ξ, δ, κ)

#             r_data = data[data .>= censoring]
#             l_data = fill(censoring, count(data .< censoring))

#             ll = sum(logcdf.(pd, l_data)) + sum(logpdf.(pd, r_data))

#             return ll
#         end
#     end

#     fobj(θ) = -loglike(θ...)

#     res = optimize(fobj, initialvalues)

#     if Optim.converged(res)
#         σ̂, ξ̂, δ̂, κ̂ = [Optim.minimizer(res)[1], Optim.minimizer(res)[2], Optim.minimizer(res)[3], Optim.minimizer(res)[4]]
#     else
#         @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#         σ̂, ξ̂, δ̂, κ̂  = [initialvalues[1], initialvalues[2], initialvalues[3], initialvalues[4]]
#     end

#     fittedmodel = EGPbetapower(σ̂, ξ̂, δ̂, κ̂)

#     return fittedmodel
# end

# function EGPbetapowerfit(data::TimeArray; initialvalues::Vector{<:Real}=Float64[], censoring::Real=0)
#     return EGPbetapowerfit(values(data), initialvalues=initialvalues, censoring=censoring)
# end


# """
#     EGPnormalfit()

# Estimate the extended GP model of Gamet and Jalbert (2020) parameters by maximum likelihood.
# """
# function EGPnormalfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])

#     if isempty(initialvalues)
#         logσ₀, ξ₀, logκ₀ = [0, 0.15, 0]  # à remplacer par une fonction (getinitialvalue)
#         initialvalues = [logσ₀, ξ₀, logκ₀]
#     end

#     function loglike(logσ::Real, ξ::Real, logκ::Real)
#         pd = EGPnormal(exp(logσ), ξ, exp(logκ))
#         ll = sum(logpdf.(pd, data))
#         return ll
#     end

#     fobj(θ) = -loglike(θ...)

#     res = optimize(fobj, initialvalues)

#     if Optim.converged(res)
#         σ̂, ξ̂, κ̂ = [exp(Optim.minimizer(res)[1]), Optim.minimizer(res)[2], exp(Optim.minimizer(res)[3])]
#     else
#         @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
#         σ̂, ξ̂, κ̂  = [exp(initialvalues[1]), initialvalues[2], exp(initialvalues[3])]
#     end

#     fittedmodel = EGPnormal(σ̂, ξ̂, κ̂)

#     return fittedmodel
# end

# function EGPnormalfit(data::TimeArray; initialvalues::Vector{<:Real}=Float64[])
#     return EGPnormalfit(values(data), initialvalues=initialvalues)
# end