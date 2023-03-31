# Application to summer precipitation recorded in Montréal, QC

In this section, we show that the EGP model based on the truncated Beta distribution can be used to model non-zero precipitation, which corresponds to exceedances above the very low threshold of 0, while maintaining the tail behavior.

```@setup rain
using Extremes, ExtendedExtremes, DataFrames, Dates, CSV, Distributions, Gadfly
```



## Data

The daily summer precipitation (May–October) recorded at the Montréal-Trudeau International Airport meteorological station (Québec, Canada) from 2000 to 2020 are investigated. This dataset can be loaded using the [`ExtendedExtremes.dataset`](@ref) provided for this tutorial.

```@example rain
data = ExtendedExtremes.dataset("pcp")
filter!(row -> row.Date>= Date(2000,1,1), data)
filter!(row -> month(row.Date) in 5:10, data)
dropmissing!(data)
first(data,5)
```

## Modeling the non-zero precipitation

The EGP parameter estimation with maximum likelihood is performed with the [`fit_mle`](@ref) function.

```@example rain
u = 0.0
y = data.pcp[data.pcp .> u] .- u;

fd = fit_mle(ExtendedGeneralizedPareto{TBeta}, y)
```

Several diagnostic plots for assessing the accuracy of the EGP model fitted to the Montréal data are can be shown with the [`diagnosticplots`](@ref) function:

```@example rain
set_default_plot_size(16cm, 16cm)
ExtendedExtremes.diagnosticplots(y, fd)
```

The diagnostic plots consist in the probability plot (upper left panel), the quantile plot (upper right panel), the density plot (lower left panel) and the return level plot (lower right panel). These plots can be displayed separately using respectively the [`probplot`](@ref), [`qqplot`](@ref), [`histplot`](@ref) and [`returnlevelplot`](@ref) functions.