"""
    dataset(name::String)::DataFrame
Load the dataset associated with `name`.
Some datasets used by Gamet & Jalbert (2022) are available using the following names:
 - `pcp`: daily rainfall accumulations in Montréal
 - `tasmax`: daily maximum temperatures recorded in Montréal
```
"""
function dataset(name::String)::DataFrame

    filename = joinpath(dirname(@__FILE__), "..", "data", string(name, ".csv"))
    if isfile(filename)
        return CSV.read(filename, DataFrame)
    end
    error("There is no dataset with the name '$name'")

end
