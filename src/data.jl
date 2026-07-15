"""
    dataset(name::String) -> DataFrame

Load the dataset associated with `name`.

The following datasets used by Gamet & Jalbert (2022) are available:

- `"pcp"`: daily rainfall accumulations in Montréal
- `"tasmax"`: daily maximum temperatures recorded in Montréal
"""
function dataset(name::String)

    filename = joinpath(@__DIR__, "..", "data", "$name.csv")

    isfile(filename) || throw(ArgumentError("There is no dataset with the name '$name'"))

    return CSV.read(filename, DataFrame)

end
