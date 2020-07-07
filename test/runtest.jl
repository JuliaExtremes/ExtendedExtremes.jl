using EGPD
using Extremes, Random, Test

# Set the seed for reproductible test results
Random.seed!(12)

@testset "EGPD.jl" begin
    include("reproducingNaveauResults.jl")
end;