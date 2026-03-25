#!/usr/bin/env julia
#
# Run benchmarks for ExtendedExtremes.
#
# Usage:
#   julia benchmark/run.jl              # run all benchmarks and print results
#   julia benchmark/run.jl --save       # run and save results to benchmark/results.json
#   julia benchmark/run.jl --compare    # run and compare against saved baseline
#

using Pkg
Pkg.activate(@__DIR__)

# Ensure BenchmarkTools is available
try
    @eval using BenchmarkTools
catch
    Pkg.add("BenchmarkTools")
    @eval using BenchmarkTools
end

include("benchmarks.jl")

const RESULTS_FILE = joinpath(@__DIR__, "results.json")

function main()
    save    = "--save"    in ARGS
    compare = "--compare" in ARGS

    println("Running ExtendedExtremes benchmarks…")
    results = run(SUITE, verbose=true)

    if save
        BenchmarkTools.save(RESULTS_FILE, median(results))
        println("\nBaseline saved to $RESULTS_FILE")
    end

    if compare && isfile(RESULTS_FILE)
        baseline = BenchmarkTools.load(RESULTS_FILE)[1]
        j = judge(median(results), baseline)
        println("\n", j)
    elseif compare
        @warn "No baseline found at $RESULTS_FILE — run with --save first."
    end

    display(median(results))
    println()
end

main()
