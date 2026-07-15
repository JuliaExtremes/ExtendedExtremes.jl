using ExtendedExtremes 
using DataFrames, Random, Test

# Set the seed for reproductible test results
Random.seed!(12)

@testset "ExtendedExtremes.jl" begin
    include("data_test.jl")
end;
