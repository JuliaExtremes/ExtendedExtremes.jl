## Quantile Plot

"""
    qqbuild()
"""
function qqbuild(d::UnivariateDistribution, x::AbstractVector)
    n = length(x)
    ind = collect(1:n)
    grid = ind ./ (n+1)

    qd = quantile.(Ref(d), grid)
    qx = sort(x, rev=false)

    return Distributions.QQPair(qd, qx)   # QQpair(x=Model, y=Empirical)
end

"""
    qqplot()

Quantile-Quantile Plot
"""
function qqplot(qq::QQPair, elements::Gadfly.ElementOrFunction...)
    return plot(x=qq.qx,  # Model
                y=qq.qy,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                #Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                #Guide.title("Quantile Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

function qqplot(fm::UnivariateDistribution, data::AbstractVector, elements::Gadfly.ElementOrFunction...)
    qq = qqbuild(fm, data)
    return plot(x=qq.qx,  # Model
                y=qq.qy,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                Guide.title("Quantile Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

## Non-stationnary Quantile Plot (À RETRAVAILLER)

struct QQPairNS{U<:AbstractVector, V<:AbstractVector}
    qx::U
    qy::V
end

"""
    qqbuild() (standardized)
"""
function qqbuild(fm::Array{<:UnivariateDistribution,1}, data::Array{<:Real,1})
    ztilde = -log.(1 .- cdf.(fm, data));  # Transformation to a standard exponential distribution
    ztildei = sort(ztilde)

    m = length(ztilde)
    x = -log.(1 .- ((1:m) / (m+1)));

    return QQPairNS(x, ztildei)
end

function qqplot(qq::QQPairNS, elements::Gadfly.ElementOrFunction...)
    return plot(x=qq.qx,  # Model
                y=qq.qy,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                #Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                #Guide.title("Residual Quantile Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

function qqplot(fm::Array{<:UnivariateDistribution,1}, data::Array{<:Real,1}, elements::Gadfly.ElementOrFunction...)
    qq = qqbuild(fm, data)
    return plot(x=qq.qx,  # Model
                y=qq.qy,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                Guide.title("Residual Quantile Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

## Probability Plot

struct PPPair{U<:AbstractVector, V<:AbstractVector}
    px::U
    py::V
end

"""
    ppbuild()
"""
function ppbuild(d::UnivariateDistribution, x::AbstractVector)

    n = length(x)
    ind = collect(1:n)
    px = ind ./ (n+1)
    pd = cdf.(Ref(d), sort(x, rev=false))

    return PPPair(pd, px)
end

"""
    probplot()

Probability plot
"""
function probplot(pp::PPPair, elements::Gadfly.ElementOrFunction...)
    return plot(x=pp.px,  # Model
                y=pp.py,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                #Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                #Guide.title("Probability Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

function probplot(fm::UnivariateDistribution, data::AbstractVector, elements::Gadfly.ElementOrFunction...)
    pp = ppbuild(fm, data)
    return plot(x=pp.px,  # Model
                y=pp.py,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                Guide.title("Probability Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

## Non-stationnary Probability Plot (À RETRAVAILLER)

struct PPPairNS{U<:AbstractVector, V<:AbstractVector}
    px::U
    py::V
end

"""
    ppbuild() (standardized)
"""
function ppbuild(fm::Array{<:UnivariateDistribution,1}, data::Array{<:Real,1})
    # TO-DO : - vérifier l'ordre des variables (x vs y)

    ztilde = -log.(1 .- cdf.(fm, data));  # Transformation to a standard exponential distribution
    ztildei = sort(ztilde)

    m = length(ztilde)
    x = (1:m) ./ (m+1);
    y = 1 .- exp.(-ztildei)

    return PPPairNS(y, x)
end

function probplot(pp::PPPairNS, elements::Gadfly.ElementOrFunction...)
    return plot(x=pp.px,  # Model
                y=pp.py,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                #Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                #Guide.title("Residual Probability Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

function probplot(fm::Array{<:UnivariateDistribution,1}, data::Array{<:Real,1}, elements::Gadfly.ElementOrFunction...)
    pp = ppbuild(fm, data)
    return plot(x=pp.px,  # Model
                y=pp.py,  # Empirical
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Guide.xlabel("Model"), Guide.ylabel("Empirical"),
                Guide.title("Residual Probability Plot"),
                Theme(discrete_highlight_color=c->nothing),
                elements...)
end

## Histogram Plot

"""
    histplot(fm::UnivariateDistribution, data::AbstractVector)

Histogram plot
"""
function histplot(fm::UnivariateDistribution, data::AbstractVector)

    #TODO: arranger max/min pour la grille

    nbin = floor(sqrt(length(data)))
    d = layer(x -> pdf(fm, x), 0, maximum(data), Geom.line, Theme(default_color="red"))
    h = layer(x = data, Geom.histogram(density=true, bincount=nbin))
    return plot(d,
                h,
                Coord.cartesian(xmin = minimum(data), xmax = maximum(data)),
                Guide.xlabel("Data"),
                Guide.ylabel("Density"),
                Guide.title("Density Plot")
                )
end

## Return Level Plot

function ecdf(y::Vector{<:Real})::Tuple{Vector{<:Real}, Vector{<:Real}}
    ys = sort(y)
    n = length(ys)
    p = collect(1:n)/(n+1)

    return ys, p
end

function returnlevelplot_data(fm::UnivariateDistribution, data::AbstractVector)
    y, p = ecdf(data)

    T = 1 ./ (1 .- p)

    n = length(y)

    q = Vector{Float64}(undef, n)

    for i in eachindex(p)
       q[i] = quantile(fm, p[i])[1]
    end

    return y, T, q

end

"""
    returnlevelplot(fm::UnivariateDistribution, data::AbstractVector)

Return level plot
"""
function returnlevelplot(fm::UnivariateDistribution, data::AbstractVector)
    y, T, q = returnlevelplot_data(fm, data)
    l1 = layer(x=T, y=q,Geom.line, Theme(default_color="red", line_style=[:dash]))
    l2 = layer(x=T, y=y, Geom.point)
    return plot(l1,
                l2,
                Scale.x_log10,
                Guide.xlabel("Return Period"),
                Guide.ylabel("Return Level"),
                Guide.title("Return Level Plot"),
                Theme(discrete_highlight_color=c->nothing))
end

## Diagnostic plots

"""
    diagnosticplots(fm::UnivariateDistribution, data::AbstractVector)

"""
function diagnosticplots(fm::UnivariateDistribution, data::AbstractVector)::Gadfly.Compose.Context
    f1 = probplot(fm, data)
    f2 = qqplot(fm, data)
    f3 = histplot(fm, data)
    f4 = returnlevelplot(fm, data)

    return gridstack([f1 f2; f3 f4])
end
