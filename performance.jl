using Pkg
pkg"activate ."

using Distributions, ExtendedExtremes, Random


pd = TBeta(.8)
y = rand(Random.seed!(1), pd, 5000)

@time  logpdf(pd, y);@time logpdf(pd, y);



pd = ExtendedGeneralizedPareto(TBeta(.8), GeneralizedPareto(1,0))
y = rand(Random.seed!(1), pd, 5000)

@time logpdf(pd, y); @time logpdf(pd, y);

@time fd = fit_mle(ExtendedGeneralizedPareto{TBeta}, y); @time fit_mle(ExtendedGeneralizedPareto{TBeta}, y);