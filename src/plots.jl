"""
    qqbuild()
"""
function qqbuild(d::UnivariateDistribution, x::AbstractVector)
    n = length(x)
    ind = collect(1:n)
    grid = ind ./ (n+1)

    qd = quantile.(Ref(d), grid)
    qx = sort(x, rev=false)

    return Distributions.QQPair(qd, qx)
end

"""
    qqbuild() (standardisé)
"""
function qqbuild(fm::Array{EGPpower{Float64},1}, data::Array{<:Real,1})
    # TO-DO : - vérifier l'ordre des variables (x vs y)

    ztilde = -log.(1 .- cdf.(fm, data));  # Transformation to a standard exponential distribution
    ztildei = sort(ztilde)

    m = length(ztilde)
    x = -log.(1 .- ((1:m) / (m+1)));

    return Distributions.QQPair(ztildei, x)
end

struct PPPair{U<:AbstractVector, V<:AbstractVector}
    px::U
    py::V
end

"""
    qqplot()
"""
function qqplot(qq::QQPair, elements::Gadfly.ElementOrFunction...)
    return plot(x=qq.qy,
                y=qq.qx,
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Theme(highlight_width=0px),
                elements...)
end


"""
    ppbuild() (standardisé)
"""
function ppbuild(fm::Array{EGPpower{Float64},1}, data::Array{<:Real,1})
    # TO-DO : - vérifier l'ordre des variables (x vs y)

    ztilde = -log.(1 .- cdf.(fm, data));  # Transformation to a standard exponential distribution
    ztildei = sort(ztilde)

    m = length(ztilde)
    x = (1:m) ./ (m+1);
    y = 1 .- exp.(-ztildei)

    return PPPair(y, x)
end

"""
    ppplot()
"""
function ppplot(pp::PPPair, elements::Gadfly.ElementOrFunction...)
    return plot(x=pp.py,
                y=pp.px,
                Geom.point,
                Geom.abline(color="red", style=:dash),
                Theme(highlight_width=0px),
                elements...)
end
