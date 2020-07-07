module EGPD

using CSV, Distributions, Extremes, Optim
import Distributions: cdf, sampler
import Random: AbstractRNG
import Base: rand

export
    # distribution types
    EGPpower
    #EGPpowermix
    #EGPbeta
    #EGPbetapower

### source files

# implementation helpers
include("utils.jl")

# distribution
include("naveau2016_type1.jl");
include("naveau2016_type2.jl");
include("naveau2016_type3.jl");
include("naveau2016_type4.jl");

end # module
