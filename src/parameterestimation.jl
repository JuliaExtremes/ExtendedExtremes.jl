function EGPpowerfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])
    inner_optimizer = NelderMead()
    lower = [0, -Inf, 0]
    upper = [Inf, Inf, Inf]

    if isempty(initialvalues)
        σ₀, ξ₀, κ₀ = [1, 0.15, 2]  # à remplacer par une fonction (getinitialvalue)
        initialvalues = [σ₀, ξ₀, κ₀]
    end

    function loglike(σ::Real, ξ::Real, κ::Real)
        pd = EGPpower(σ, ξ, κ)
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, lower, upper, initialvalues, Fminbox(inner_optimizer))

    if Optim.converged(res)
        σ̂, ξ̂, κ̂ = [Optim.minimizer(res)...]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, κ̂  = [initialvalues...]
    end

    fittedmodel = EGPpower(σ̂, ξ̂, κ̂)

    return fittedmodel
end

function EGPpowermixfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])
    inner_optimizer = NelderMead()
    lower = [0, -Inf, 0, 0, 0]
    upper = [Inf, Inf, Inf, Inf, 1]

    if isempty(initialvalues)
        σ₀, ξ₀, κ₁₀, κ₂₀, p₀ = [1, .15, 1, 2, .4]  # à remplacer par une fonction (getinitialvalue)
        initialvalues = [σ₀, ξ₀, κ₁₀, κ₂₀, p₀]
    end

    function loglike(σ::Real, ξ::Real, κ₁::Real, κ₂::Real, p::Real)
        pd = EGPpowermix(σ, ξ, κ₁, κ₂, p)
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, lower, upper, initialvalues, Fminbox(inner_optimizer))

    if Optim.converged(res)
        σ̂, ξ̂, κ̂₁, κ̂₂, p̂ = [Optim.minimizer(res)...]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, κ̂₁, κ̂₂, p̂  = [initialvalues...]
    end

    fittedmodel = EGPpowermix(σ̂, ξ̂, κ̂₁, κ̂₂, p̂)

    return fittedmodel
end

function EGPbetafit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])
    inner_optimizer = NelderMead()
    lower = [0, -Inf, 0]
    upper = [Inf, Inf, Inf]

    if isempty(initialvalues)
        σ₀, ξ₀, δ₀ = [1, 0.15, 5]  # à remplacer par une fonction (getinitialvalue)
        initialvalues = [σ₀, ξ₀, δ₀]
    end

    function loglike(σ::Real, ξ::Real, δ::Real)
        pd = EGPbeta(σ, ξ, δ)
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, lower, upper, initialvalues, Fminbox(inner_optimizer))

    if Optim.converged(res)
        σ̂, ξ̂, δ̂ = [Optim.minimizer(res)...]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, δ̂  = [initialvalues...]
    end

    fittedmodel = EGPbeta(σ̂, ξ̂, δ̂)

    return fittedmodel
end

function EGPbetapowerfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])
    inner_optimizer = NelderMead()
    lower = [0, -Inf, 0, 0]
    upper = [Inf, Inf, Inf, Inf]

    if isempty(initialvalues)
        σ₀, ξ₀, δ₀, κ₀ = [1, 0.15, 2, 2]  # à remplacer par une fonction (getinitialvalue)
        initialvalues = [σ₀, ξ₀, δ₀, κ₀]
    end

    function loglike(σ::Real, ξ::Real, δ::Real, κ::Real)
        pd = EGPbetapower(σ, ξ, δ, κ)
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, lower, upper, initialvalues, Fminbox(inner_optimizer))

    if Optim.converged(res)
        σ̂, ξ̂, δ̂, κ̂ = [Optim.minimizer(res)...]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, δ̂, κ̂  = [initialvalues...]
    end

    fittedmodel = EGPbetapower(σ̂, ξ̂, δ̂, κ̂)

    return fittedmodel
end

function EGPnormalfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])

    if isempty(initialvalues)
        logσ₀, ξ₀, logκ₀ = [0, 0.15, 0]  # à remplacer par une fonction (getinitialvalue)
        initialvalues = [logσ₀, ξ₀, logκ₀]
    end

    function loglike(logσ::Real, ξ::Real, logκ::Real)
        pd = EGPnormal(exp(logσ), ξ, exp(logκ))
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        σ̂, ξ̂, κ̂ = [exp(Optim.minimizer(res)[1]), Optim.minimizer(res)[2], exp(Optim.minimizer(res)[3])]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, κ̂  = [exp(initialvalues[1]), initialvalues[2], exp(initialvalues[3])]
    end

    fittedmodel = EGPnormal(σ̂, ξ̂, κ̂)

    return fittedmodel
end
