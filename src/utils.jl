# macro for argument checking (taken from Distributions.jl)

macro check_args(D, cond)
    quote
        if !($(esc(cond)))
            throw(ArgumentError(string(
                $(string(D)), ": the condition ", $(string(cond)), " is not satisfied.")))
        end
    end
end

"""
    BIC(fm::ContinuousUnivariateDistribution, data::Array{<:Real,1}; censoring::Real=0)

Calculates the Bayesian information criterion (BIC) of the specified model.
"""
function BIC(fm::ContinuousUnivariateDistribution, data::Array{<:Real,1}; censoring::Real=0)
    n = size(data,1)
    k = size(params(fm),1)

    r_data = data[data .>= censoring]
    l_data = fill(censoring, count(data .< censoring))

    ll = sum(logcdf.(fm, l_data)) + sum(logpdf.(fm, r_data))

    return ll - k/2*log(n)
end

function BIC(fm::ContinuousUnivariateDistribution, data::TimeArray; censoring::Real=0)
    return BIC(fm, values(data), censoring=censoring)
end

function window(ts::TimeArray, year_from::Int64, year_to::Int64)
    return to(from(ts, DateTime(year_from, 01, 01)), DateTime(year_to-1, 12, 31))
end
