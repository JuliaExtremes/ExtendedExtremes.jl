@testset "EGPbeta" begin
    @test_throws ArgumentError EGPbeta(-1, .2, 2)
    @test_throws ArgumentError EGPbeta(-1, .2, -2)

    σ = 1
    ξ = .3
    δ = 2
    model = EGPpower(σ, ξ, δ)

    @test minimum(model) == 0
    @test maximum(model) == Inf

    @test insupport(model, -1) == false
    @test insupport(model, 0) == true
    @test insupport(model, 1000) == true

    @test params(model) == (σ, ξ, δ)

    @test partype(model) == Float64
end
