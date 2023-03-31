module ExtendedExtremes

using Distributions
using Optim

import Distributions.@check_args
import Distributions.params
import Distributions.scale
import Distributions.shape
import Distributions.minimum
import Distributions.maximum
import Distributions.insupport
import Distributions.cdf
import Distributions.logpdf
import Distributions.pdf
import Distributions.quantile
import Distributions.rand
import Distributions.sampler
import Distributions.fit_mle
import Distributions.fit
import Random: AbstractRNG

# distributions
include("distributions/TruncatedBeta.jl");
include("distributions/TruncatedNormal.jl");
include("distributions/Power.jl");
include("distributions/ExtendedGeneralizedPareto.jl");
include("parameterestimation.jl");
include("plots.jl")

export
    # distribution types
    TBeta,
    TNormal,
    Power,
    ExtendedGeneralizedPareto,

    # methods
    EGPtype,
    params,      # get the tuple of parameters
    partype,
    scale,
    tailindex,
    shape,
    minimum,
    maximum,
    insupport,   # predicate, is x in the support of the distribution?
    getdistribution,
    #logcdf,      # cdf returning log-probability
    cdf,         # cumulative distribution function
    logpdf,      # log probability density
    pdf,         # probability density function
    quantile,    # inverse of cdf (defined for p in (0,1))
    rand,
    #sampler,      # create a Sampler object for efficient samples

    # distribution fitting
    fit_mle

end # module
