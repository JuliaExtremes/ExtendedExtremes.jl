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
    BIC(fm::ContinuousUnivariateDistribution, data::Array{<:Real,1})

Calculates the Bayesian information criterion (BIC) of the specified model.
"""
function BIC(fm::ContinuousUnivariateDistribution, data::Array{<:Real,1})
    n = size(data,1)
    k = size(params(fm),1)

    return sum(logpdf.(fm, data)) - k/2*log(n)
end
