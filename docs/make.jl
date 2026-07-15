using Documenter
using ExtendedExtremes
using ExtremePlots

using Distributions
using Random

using Gadfly
import Pango_jll
import Cairo
import Fontconfig

makedocs(modules = [ExtendedExtremes],
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
