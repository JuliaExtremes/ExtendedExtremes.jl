module EGPD

using CSV, Distributions, Extremes, Optim
import Distributions: cdf, sampler
import Random: AbstractRNG
import Base: maximum, minimum, rand

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
    params,      # get the tuple of parameters
    partype,
    pdf,         # probability density function
    quantile,    # inverse of cdf (defined for p in (0,1))
    rand,
    sampler      # create a Sampler object for efficient samples


### source files

# implementation helpers
include("utils.jl")

# distribution
include("naveau2016_type1.jl");
include("naveau2016_type2.jl");
include("naveau2016_type3.jl");
include("naveau2016_type4.jl");

end # module
