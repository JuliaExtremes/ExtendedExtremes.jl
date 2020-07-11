module EGPD

using CSV, Distributions, Extremes, Optim
import Distributions: cdf, insupport, logcdf, logpdf, params, partype, pdf, quantile, sampler
import Random: AbstractRNG
import Base: maximum, minimum, rand

include("parameterestimation.jl")
include("utils.jl")

# distributions
include("distributions/naveau2016_type1.jl");
include("distributions/naveau2016_type2.jl");
include("distributions/naveau2016_type3.jl");
include("distributions/naveau2016_type4.jl");

export
    # distribution types
    EGPpower,
    EGPpowermix,
    EGPbeta,
    EGPbetapower,

    # methods
    cdf,         # cumulative distribution function
    insupport,   # predicate, is x in the support of the distribution?
    logpdf,      # log probability density
    logcdf,      # cdf returning log-probability
    maximum,
    minimum,
    params,      # get the tuple of parameters
    partype,
    pdf,         # probability density function
    quantile,    # inverse of cdf (defined for p in (0,1))
    rand,
    sampler,      # create a Sampler object for efficient samples

    EGPpowerfit,
    EGPpowermixfit,
    EGPbetafit,
    EGPbetapowerfit

end # module
