using CSV, DataFrames, Dates, Distributions, ExtendedExtremes, Extremes, QuantileMatching


member = "kda"
filename = "/Users/jalbert/Library/CloudStorage/Dropbox/Files/Data/Simulations/CRCM/climex/ClimEx_mtl/$member.csv"

df_sim = CSV.read(filename, DataFrame)
rename!(df_sim, :Precipitation => :RAW)
filter!(row ->2000 ≤ year(row.Date) ≤ 2020, df_sim)
filter!(row -> 5 ≤ month(row.Date) ≤ 10, df_sim)
first(df_sim, 5)

p = .41

z = collect(skipmissing(df_sim.RAW))
u = wet_threshold(z, p)
x = censor(z, u)

df_results = DataFrame(Date = df_sim.Date, RAW = x)
filter!(row -> row.RAW>0, df_results)
first(df_results, 5)

x⁺ = df_results.RAW;   

fₐ = fit_mle(ExtendedGeneralizedPareto{TBeta}, x⁺, leftcensoring = 0.)

QuantileMatching.qqplot(x⁺, fₐ)




b = 1/2
V(a, α) = LocationScale(-a/(b-a), 1/(b-a), Truncated(Beta(α, α), a, b))

function fobj(θ::AbstractVector{<:Real})

    a, α, σ, ξ = θ

    if (0 < a < 1/2) & (α>0)

        pd = ExtendedGeneralizedPareto(V(a, α), GeneralizedPareto(σ, ξ))

        return -sum(logpdf.(pd, x⁺))

    else

        return Inf
    
    end

end

res = optimize(fobj, [1/8, 1, 1, 0])


â, κ̂, σ̂, ξ̂ = Optim.minimizer(res)

pd = ExtendedGeneralizedPareto(V(â,κ̂), GeneralizedPareto(σ̂, ξ̂))

QuantileMatching.qqplot(x⁺, pd)



# u = quantile(x⁺, .95)

# z = x⁺[x⁺ .> u] .-u

# pd = Extremes.getdistribution(gpfit(z))[]

# QuantileMatching.qqplot(z, pd)
