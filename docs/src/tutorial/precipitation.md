# Application to summer precipitation recorded in Montréal, QC

In this section, we show that the EGP model based on the truncated Beta distribution can be used to model non-zero precipitation, which corresponds to exceedances above the very low threshold of 0, while maintaining the tail behavior.

```@setup rain
using Extremes, ExtendedExtremes, DataFrames, Dates, CSV, Distributions, Gadfly
```



## Data

```@example rain
data = ExtendedExtremes.dataset("pcp")
filter!(row -> row.Date>= Date(2000,1,1), data)
filter!(row -> month(row.Date) in 5:10, data)
dropmissing!(data)
first(data,5)
```


## Tail analysis

```@example rain
u = 30
z = data.pcp[data.pcp.>u] .- u

fm = gpfit(z)
```

## Modeling the non-zero precipitation

```@example rain
u = 0.0
y = data.pcp[data.pcp .> u] .- u;

fd = fit_mle(ExtendedGeneralizedPareto{TBeta}, y)

set_default_plot_size(16cm, 16cm)
ExtendedExtremes.diagnosticplots(y, fd)
```