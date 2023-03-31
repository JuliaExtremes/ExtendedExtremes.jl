"""
    histplot(y::Vector{Float64}, pd::Distribution)

Histogram plot
"""
function histplot(y::Vector{Float64}, pd::Distribution)
    
    n = length(y)
    
    lbound = 0
    ubound = quantile(pd, 1-1/n/100)
    
    return histplot(y, pd, lbound, ubound)
    
end

"""
    histplot(y::Vector{Float64}, pd::Distribution, lbound::Real, ubound::Real)

Histogram plot
"""
function histplot(y::Vector{Float64}, pd::Distribution, lbound::Real, ubound::Real)
    
    n = length(y)
    nbin = Int64(ceil(sqrt(n)))
    
    x = range(lbound, stop=ubound, length=100)
    h = layer(x = y, Geom.histogram(density=true, bincount=nbin), alpha=[1.0])
    d = layer(x -> pdf(pd, x), lbound , ubound, linestyle=[:solid], Theme(default_color="red"))
    
    return plot(d, h, Coord.cartesian(xmin = lbound, xmax = ubound),
    Guide.xlabel("Data"), Guide.ylabel("Density"), Guide.title("Density plot"),
    Theme(discrete_highlight_color=c->nothing))
    
end

"""
    qqplot(y::Vector{Float64}, pd::Distribution)

Quantile-Quantile plot
"""
function qqplot(y::Vector{Float64}, pd::Distribution)

    n = length(y)
    q = sort(y)

    p = (1:n) ./ (n+1)

    q̂ = quantile.(pd, p);
    
    l1 = layer(y = q, x = q̂, Geom.point)
    l2 = layer(y = q[[1, end]], x = q[[1, end]], Geom.line, Theme(default_color="red", line_style=[:dash]))

    return plot(l2, l1,
        Guide.xlabel("Model"), Guide.ylabel("Empirical"), Guide.title("Quantile Plot"),
        Theme(discrete_highlight_color=c->nothing))
end


"""
    returnlevelplot(y::Vector{Float64}, pd::Distribution)

Return level plot
"""
function returnlevelplot(y::Vector{Float64}, pd::Distribution)

    n = length(y)
    q = sort(y)

    p = (1:n) ./ (n+1)

    t = 1 ./ (1 .- p);
    
    l1 = layer(x = t, y = q, Geom.point)
    l2 = layer(x = t, y = quantile.(pd, p), Geom.line, Theme(default_color="red", line_style=[:dash]))

    return plot(l2, l1, Scale.x_log10, Guide.xlabel("Return Period"), Guide.ylabel("Return Level"),
        Guide.title("Return Level Plot"), Theme(discrete_highlight_color=c->nothing))
        
end

"""
    probplot(y::Vector{Float64}, pd::Distribution)

Probability plot
"""
function probplot(y::Vector{Float64}, pd::Distribution)
    
    n = length(y)
    p = (1:n) ./ (n+1)
    
    q = sort(y)
    p̂ = cdf.(pd, q)
    
    l1 = layer(y = p, x = p̂, Geom.point)
    l2 = layer(y = p[[1, end]], x = p[[1, end]], Geom.line, Theme(default_color="red", line_style=[:dash]))

    return plot(l2, l1,
        Guide.xlabel("Model"), Guide.ylabel("Empirical"), Guide.title("Probability Plot"),
        Theme(discrete_highlight_color=c->nothing))
end

"""
    diagnosticplots(fm::fittedEVA)

Diagnostic plots
"""
function diagnosticplots(y::Vector{Float64}, pd::Distribution)::Gadfly.Compose.Context

	f1 = probplot(y, pd)
	f2 = qqplot(y, pd)
	f3 = histplot(y, pd)
    f4 = returnlevelplot(y, pd)

    return gridstack([f1 f2; f3 f4])
end
