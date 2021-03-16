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
