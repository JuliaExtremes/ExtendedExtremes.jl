using Documenter, Extremes, ExtendedExtremes, Distributions, Gadfly,  Random, Cairo, Fontconfig

CI = get(ENV, "CI", nothing) == "true"

makedocs(
    sitename = "ExtendedExtremes.jl",
    format = Documenter.HTML(
    prettyurls = CI,
    ),
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

if CI
    deploydocs(
    repo   = "github.com/JuliaExtremes/ExtendedExtremes.jl.git",
    devbranch = "dev",
    versions = ["stable" => "v^", "v#.#", "master"],
    push_preview = false,
    target = "build"
    )
end