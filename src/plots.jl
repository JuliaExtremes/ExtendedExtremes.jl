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
