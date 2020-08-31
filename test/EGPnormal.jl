@testset "EGPnormal" begin
    @test_throws ArgumentError EGPnormal(-1, .2, 2)
    @test_throws ArgumentError EGPnormal(-1, .2, -2)

    σ = 1
    ξ = -0.1
    κ = 2
    model = EGPnormal(σ, ξ, κ)
    @test minimum(model) == 0
    @test maximum(model) == -σ/ξ

    σ = 1
    ξ = .1
    κ = 2
    model = EGPnormal(σ, ξ, κ)

    @test maximum(model) == Inf

    @test insupport(model, -1) == false
    @test insupport(model, 0) == true
    @test insupport(model, 1000) == true

    @test params(model) == (σ, ξ, κ)

    @test partype(model) == Float64

    @test_throws AssertionError quantile(model, 0)
    @test_throws AssertionError quantile(model, 1)
end
