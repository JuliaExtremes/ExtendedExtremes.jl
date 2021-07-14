# Skewness
"""
    index_skew(ts::TimeArray)

Returns the skewness of ts.
"""
function index_skew(ts::TimeArray)
    x = values(ts)
    μ = mean(x)
    σ = std(x)

    z = (x .- μ) ./ σ

    return mean(z .^ 3)
end


"""
    index_Rx(ts::TimeArray, x::Real)

Returns the relative frequency of days with precipitation above x mm.
"""
function index_Rx(ts::TimeArray, x::Real)
    return sum(values(ts .> x)) / length(ts)
end

# R10
"""
    index_R10(ts::TimeArray)

Returns the relative frequency of days with precipitation above 10 mm.
"""
function index_R10(ts::TimeArray)
    return index_Rx(ts, 10)
end


"""
    index_pXwet(ts::TimeArray, p::Real; wet_thresh::Real=1.0)

Returns the p percentile of wet days considering days with precipitation above wet_thresh (mm).
"""
function index_pXwet(ts::TimeArray, p::Real; wet_thresh::Real=1.0)
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    x = values(ts)
    x = x[x .>= wet_thresh]

    return quantile(x, p)
end

# P98Wet
"""
    index_p98wet(ts::TimeArray)

Returns the 98% percentile of wet days considering days with precipitation above 1 mm.
"""
function index_p98wet(ts::TimeArray)
    return index_pXwet(ts, 0.98)
end

# P90Wet
"""
    index_p90wet(ts::TimeArray)

Returns the 90% percentile of wet days considering days with precipitation above 1 mm.
"""
function index_p90wet(ts::TimeArray)
    return index_pXwet(ts, 0.9)
end


"""
    index_pXwetamount(ts::TimeArray, p::Real; wet_thresh=1.0)

Returns the precipitation amount from days above p percentile.
"""
function index_pXwetamount(ts::TimeArray, p::Real; wet_thresh=1.0)
    @assert zero(p)<p<one(p) "the quantile level should be between 0 and 1."

    thresh = index_pXwet(ts, p, wet_thresh=wet_thresh)

    x = values(ts)
    x = x[x .> thresh]

    return sum(x) / length(x)
end

# P98WetAmount
"""
    index_p98wetamount(ts::TimeArray)

Returns the precipitation amount from days above 98% percentile.
"""
function index_p98wetamount(ts::TimeArray)
    return index_pXwetamount(ts, 0.98)
end


"""
    index_RVxmax(ts::TimeArray, p::Real)

Returns the p-years return value of heavy precipitation.
"""
function index_RVxmax(ts::TimeArray, p::Real)

    x = collapse(ts,year,maximum)

    fm = Extremes.gevfit(values(x))
    r = Extremes.returnlevel(fm, p)

    return r.value[1]
end

# RV20_max
"""
    index_RV20max(ts::TimeArray)

Returns the 20-years return value of heavy precipitation.
"""
function index_RV20max(ts::TimeArray)
    return index_RVxmax(ts, 20.0)
end

# SDII
"""
    index_SDII(ts::TimeArray)

Returns the precipitation amount in wet days divided by the number of wet days.
"""
function index_SDII(ts::TimeArray; thresh::Real=1.0)

    x = values(ts)
    idx = x .> thresh

    return sum(x[idx]) / length(x[idx])
end
