using Documenter, ExtendedExtremes, Distributions, Random

makedocs(modules = [ExtendedExtremes, Distributions, Random],
        doctest = false,
        sitename="ExtendedExtremes.jl",
        pages = [
        "index.md",
        "starting.md",
        "distributions.md",
        "fit.md",
        "diagnostics.md",
        "contributing.md",
        #"Tutorial" =>["Getting started" => "tutorial/index.md"],
        #"contributing.md",
        "functions.md"
        ]
)

deploydocs(
        repo = "github.com/houton199/ExtendedExtremes.jl.git",
)
