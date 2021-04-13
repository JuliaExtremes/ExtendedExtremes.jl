module EGPD

using CSV, Distributions, Extremes, Optim
using Gadfly
import Distributions: cdf, insupport, logcdf, logpdf, params, partype, pdf, quantile, sampler
#import Distributions: qqbuild
import Random: AbstractRNG
import Base: maximum, minimum, rand

include("parameterestimation.jl")
include("plots.jl")
include("utils.jl")

# distributions
include("distributions/naveau2016_type1.jl");
include("distributions/naveau2016_type1_nonstat.jl");
include("distributions/naveau2016_type2.jl");
include("distributions/naveau2016_type3.jl");
include("distributions/naveau2016_type4.jl");
include("distributions/EGPnormal.jl");

# bias correction
include("biascorrect.jl");

export
    # distribution types
    EGPpower,
    EGPpowermix,
    EGPbeta,
    EGPbetapower,
    EGPnormal,

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

    #load,

    BIC,

    qqbuild,
    qqplot,
    Guide,

    EGPpowerfit,
    EGPpowermixfit,
    EGPbetafit,
    EGPbetapowerfit,
    EGPnormalfit,

    T_transform,
    #invT_transform
    bias_correct

end # module
