function EGPpowerfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])

    if isempty(initialvalues)
        logσ₀, ξ₀, logκ₀ = [0, 0.15, 0]
        initialvalues = [logσ₀, ξ₀, logκ₀]
    end

    function loglike(logσ::Real, ξ::Real, logκ::Real)
        pd = EGPpower(exp(logσ), ξ, exp(logκ))
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

    fittedmodel = EGPpower(σ̂, ξ̂, κ̂)

    return fittedmodel
end

function EGPbetafit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])

    if isempty(initialvalues)
        logσ₀, ξ₀, logδ₀ = [0, 0.15, 0]
        initialvalues = [logσ₀, ξ₀, logδ₀]
    end

    function loglike(logσ::Real, ξ::Real, logδ::Real)
        pd = EGPbeta(exp(logσ), ξ, exp(logδ))
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        σ̂, ξ̂, δ̂ = [exp(Optim.minimizer(res)[1]), Optim.minimizer(res)[2], exp(Optim.minimizer(res)[3])]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, δ̂  = [exp(initialvalues[1]), initialvalues[2], exp(initialvalues[3])]
    end

    fittedmodel = EGPbeta(σ̂, ξ̂, δ̂)

    return fittedmodel
end

function EGPbetapowerfit(data::Array{<:Real,1}; initialvalues::Vector{<:Real}=Float64[])

    if isempty(initialvalues)
        logσ₀, ξ₀, logδ₀, logκ₀ = [0, 0.15, 0, 0]
        initialvalues = [logσ₀, ξ₀, logδ₀, logκ₀]
    end

    function loglike(logσ::Real, ξ::Real, logδ::Real, logκ::Real)
        pd = EGPbetapower(exp(logσ), ξ, exp(logδ), exp(logκ))
        ll = sum(logpdf.(pd, data))
        return ll
    end

    fobj(θ) = -loglike(θ...)

    res = optimize(fobj, initialvalues)

    if Optim.converged(res)
        σ̂, ξ̂, δ̂, κ̂ = [exp(Optim.minimizer(res)[1]), Optim.minimizer(res)[2], exp(Optim.minimizer(res)[3]), exp(Optim.minimizer(res)[4])]
    else
        @warn "The maximum likelihood algorithm did not find a solution. Maybe try with different initial values or with another method. The returned values are the initial values."
        σ̂, ξ̂, δ̂, κ̂  = [exp(initialvalues[1]), initialvalues[2], exp(initialvalues[3]), exp(initialvalues[4])]
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
