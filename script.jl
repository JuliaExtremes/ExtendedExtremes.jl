using Distributions, ExtendedExtremes

pd = ExtendedGeneralizedPareto(TBeta(.5), GeneralizedPareto(1, 0))

y = rand(pd, 10000)

fit_mle(ExtendedGeneralizedPareto{TBeta}, y, [.5, 1, 0], leftcensoring = .01)

leftcensoring = 0.

y⁺ = filter( v -> v > leftcensoring, y)

n⁻ = count( y .< leftcensoring)

sum(logpdf.(pd, y⁺)) + n⁻ * logcdf(pd, leftcensoring)


