using BenchmarkTools
using ExtendedExtremes
using Distributions
using Random

Random.seed!(42)

const SUITE = BenchmarkGroup()

# ============================================================================
# Distribution construction
# ============================================================================
SUITE["construction"] = BenchmarkGroup()
SUITE["construction"]["Power"]          = @benchmarkable Power(2.0)
SUITE["construction"]["TNormal"]        = @benchmarkable TNormal(3.0)
SUITE["construction"]["TBeta"]          = @benchmarkable TBeta(2.0)
SUITE["construction"]["EGP_Power"]      = @benchmarkable ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1))
SUITE["construction"]["EGP_TNormal"]    = @benchmarkable ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1))
SUITE["construction"]["EGP_TBeta"]      = @benchmarkable ExtendedGeneralizedPareto(TBeta(2.0), GeneralizedPareto(1.0, 0.1))

# ============================================================================
# PDF evaluation
# ============================================================================
SUITE["pdf"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["pdf"]["EGP_Power_scalar"]    = @benchmarkable pdf($pd, 1.5)
    SUITE["pdf"]["EGP_Power_vector"]    = @benchmarkable pdf.($pd, $x)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["pdf"]["EGP_TNormal_scalar"]  = @benchmarkable pdf($pd, 1.5)
    SUITE["pdf"]["EGP_TNormal_vector"]  = @benchmarkable pdf.($pd, $x)
end

let pd = ExtendedGeneralizedPareto(TBeta(2.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["pdf"]["EGP_TBeta_scalar"]    = @benchmarkable pdf($pd, 1.5)
    SUITE["pdf"]["EGP_TBeta_vector"]    = @benchmarkable pdf.($pd, $x)
end

# ============================================================================
# logpdf evaluation
# ============================================================================
SUITE["logpdf"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["logpdf"]["EGP_Power_scalar"]  = @benchmarkable logpdf($pd, 1.5)
    SUITE["logpdf"]["EGP_Power_vector"]  = @benchmarkable logpdf.($pd, $x)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["logpdf"]["EGP_TNormal_scalar"] = @benchmarkable logpdf($pd, 1.5)
    SUITE["logpdf"]["EGP_TNormal_vector"] = @benchmarkable logpdf.($pd, $x)
end

# ============================================================================
# CDF evaluation
# ============================================================================
SUITE["cdf"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["cdf"]["EGP_Power_scalar"]    = @benchmarkable cdf($pd, 1.5)
    SUITE["cdf"]["EGP_Power_vector"]    = @benchmarkable cdf.($pd, $x)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1)),
    x = rand(pd, 1000)
    SUITE["cdf"]["EGP_TNormal_scalar"]  = @benchmarkable cdf($pd, 1.5)
    SUITE["cdf"]["EGP_TNormal_vector"]  = @benchmarkable cdf.($pd, $x)
end

# ============================================================================
# Quantile evaluation
# ============================================================================
SUITE["quantile"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    p = collect(range(0.01, 0.99, length=100))
    SUITE["quantile"]["EGP_Power_scalar"]  = @benchmarkable quantile($pd, 0.95)
    SUITE["quantile"]["EGP_Power_vector"]  = @benchmarkable quantile.($pd, $p)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1)),
    p = collect(range(0.01, 0.99, length=100))
    SUITE["quantile"]["EGP_TNormal_scalar"] = @benchmarkable quantile($pd, 0.95)
    SUITE["quantile"]["EGP_TNormal_vector"] = @benchmarkable quantile.($pd, $p)
end

# ============================================================================
# Random sampling
# ============================================================================
SUITE["rand"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1))
    SUITE["rand"]["EGP_Power_single"]   = @benchmarkable rand($pd)
    SUITE["rand"]["EGP_Power_1000"]     = @benchmarkable rand($pd, 1000)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1))
    SUITE["rand"]["EGP_TNormal_single"] = @benchmarkable rand($pd)
    SUITE["rand"]["EGP_TNormal_1000"]   = @benchmarkable rand($pd, 1000)
end

let pd = ExtendedGeneralizedPareto(TBeta(2.0), GeneralizedPareto(1.0, 0.1))
    SUITE["rand"]["EGP_TBeta_single"]   = @benchmarkable rand($pd)
    SUITE["rand"]["EGP_TBeta_1000"]     = @benchmarkable rand($pd, 1000)
end

# ============================================================================
# Parameter estimation (fit_mle)
# ============================================================================
SUITE["fit_mle"] = BenchmarkGroup()

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    data = rand(pd, 500)
    SUITE["fit_mle"]["EGP_Power_n500"] = @benchmarkable fit_mle(ExtendedGeneralizedPareto{Power}, $data)
end

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    data = rand(pd, 2000)
    SUITE["fit_mle"]["EGP_Power_n2000"] = @benchmarkable fit_mle(ExtendedGeneralizedPareto{Power}, $data)
end

let pd = ExtendedGeneralizedPareto(TNormal(3.0), GeneralizedPareto(1.0, 0.1)),
    data = rand(pd, 500)
    SUITE["fit_mle"]["EGP_TNormal_n500"] = @benchmarkable fit_mle(ExtendedGeneralizedPareto{TNormal}, $data)
end

let pd = ExtendedGeneralizedPareto(Power(2.0), GeneralizedPareto(1.0, 0.1)),
    data = rand(pd, 500)
    SUITE["fit_mle"]["EGP_Power_censored_n500"] = @benchmarkable fit_mle(ExtendedGeneralizedPareto{Power}, $data; leftcensoring=0.5)
end

# ============================================================================
# Data loading
# ============================================================================
SUITE["data"] = BenchmarkGroup()
SUITE["data"]["load_pcp"]    = @benchmarkable ExtendedExtremes.dataset("pcp")
SUITE["data"]["load_tasmax"] = @benchmarkable ExtendedExtremes.dataset("tasmax")

# ============================================================================
# End-to-end workflow: load data → fit → quantile
# ============================================================================
SUITE["workflow"] = BenchmarkGroup()

let pcp = ExtendedExtremes.dataset("pcp"),
    y = collect(skipmissing(filter(x -> !ismissing(x) && x > 0, pcp[!, end])))
    SUITE["workflow"]["pcp_fit_and_quantile"] = @benchmarkable begin
        fd = fit_mle(ExtendedGeneralizedPareto{Power}, $y)
        quantile(fd, 0.99)
    end
end
