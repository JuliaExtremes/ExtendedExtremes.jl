module EGPD

using CSV, Extremes, Optim
import Distributions

include("naveau2016_type1.jl");
include("naveau2016_type2.jl");
include("naveau2016_type3.jl");
include("naveau2016_type4.jl");

end # module
