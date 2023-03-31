# Extended generalized Pareto (EGP) distributions

The extended generalized Pareto families proposed in [Gamet & Jalbert (2022)](https://doi.org/10.1002/env.2744) are implemented and ready to be used just like any [`UnivariateDistribution`](https://juliastats.org/Distributions.jl/stable/univariate/) from the *[Distributions.jl](https://github.com/JuliaStats/Distributions.jl)* package.


## EGP model based on the power function


```@setup power
using ExtendedExtremes, Distributions, Gadfly
```

```@example power
l1 = layer(x -> pdf(ExtendedGeneralizedPareto(Power(0.8), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 0.8"], linestyle=[:dash])
l2 = layer(x -> pdf(ExtendedGeneralizedPareto(Power(1.0), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 1.0"])
l3 = layer(x -> pdf(ExtendedGeneralizedPareto(Power(1.5), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 1.5"], linestyle=[:dot])

xticks = [0.00, 0.50, 1.00, 1.50, 2.00]
yticks = [0.00, 0.50, 1.00, 1.50, 2.00]

set_default_plot_size(13cm, 10cm)
plot(l1, l2, l3,
    Guide.colorkey(title="", labels=["κ = 0.8","κ = 1.0", "κ = 1.5"]),
    Guide.xticks(ticks=xticks), Guide.yticks(ticks=yticks),
    Guide.ylabel(""), Guide.xlabel(""), Guide.title("EGP-Power")
)
```


## EGP model based on the truncated normal distribution

```@setup normal
using ExtendedExtremes, Distributions, Gadfly
```

```@example normal
l1 = layer(x -> pdf(ExtendedGeneralizedPareto(TNormal(.00001), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 0"], linestyle=[:dash])
l2 = layer(x -> pdf(ExtendedGeneralizedPareto(TNormal(1), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 1"])
l3 = layer(x -> pdf(ExtendedGeneralizedPareto(TNormal(5), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 5"], linestyle=[:dot])

xticks = [0.00, 0.50, 1.00, 1.50, 2.00]
yticks = [0.20, 0.40, 0.60, 0.80, 1.00]

set_default_plot_size(13cm, 10cm)
plot(l1, l2, l3,
    Guide.colorkey(title="", labels=["κ = 0","κ = 1", "κ = 5"]),
    Guide.xticks(ticks=xticks), Guide.yticks(ticks=yticks),
    Guide.ylabel(""), Guide.xlabel(""), Guide.title("EGP-TNormal")
)
```


## EGP model based on the truncated Beta distribution

```@setup beta
using ExtendedExtremes, Distributions, Gadfly
```

```@example beta
l1 = layer(x -> pdf(ExtendedGeneralizedPareto(TBeta(0.8), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 0.8"], linestyle=[:dash])
l2 = layer(x -> pdf(ExtendedGeneralizedPareto(TBeta(1.0), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 1.0"])
l3 = layer(x -> pdf(ExtendedGeneralizedPareto(TBeta(1.5), GeneralizedPareto(1,0)), x), 0 , 2, color=["κ = 1.5"], linestyle=[:dot])

xticks = [0.00, 0.50, 1.00, 1.50, 2.00]
yticks = [0.25, 0.50, 0.75, 1.00, 1.25]

set_default_plot_size(13cm, 10cm)
plot(l1, l2, l3,
    Guide.colorkey(title="", labels=["κ = 0.8","κ = 1.0", "κ = 1.5"]),
    Guide.xticks(ticks=xticks), Guide.yticks(ticks=yticks),
    Guide.ylabel(""), Guide.xlabel(""), Guide.title("EGP-TBeta")
)
```
