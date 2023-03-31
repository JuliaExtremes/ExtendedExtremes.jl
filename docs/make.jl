using Documenter, ExtendedExtremes, Distributions, Random, Cairo, Fontconfig

makedocs(#modules = [ExtendedExtremes, Distributions, Random, Cairo, Fontconfig],
        doctest = false,
        sitename="ExtendedExtremes.jl",
        pages = [
        "index.md",
		"Tutorial" =>["Getting started" => "tutorial/index.md",
					"Extended GP distributions" => "tutorial/distributions.md",
					"Application: Precipitation" => "tutorial/precipitation.md"],
					#"Application: Temperatures" => "tutorial/temperatures.md"],
        "contributing.md",
        "functions.md"
        ]
)

deploydocs(
        repo = "github.com/JuliaExtremes/ExtendedExtremes.jl.git",
)
