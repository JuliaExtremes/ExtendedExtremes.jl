@testset "Reproducing Naveau Results" begin
    @testset "Model (i)" begin

        data = Extremes.dataset("rain")

        # Fit of the first model by maximum likelihood
        fm = EGPpowerfit(data.Rainfall[data.Rainfall .> 0])

        # Parameter estimates
        θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]]

        # Parameter estimates with Naveau (R)
        θ = [4.56372682168789; 0.223166954694046; 1.19212081334056]

        @test θ̂ ≈ θ rtol = 0.01
    end

    @testset "Model (iii)" begin

        data = Extremes.dataset("rain")

        # Fit of the first model by maximum likelihood
        fm = EGPbetafit(data.Rainfall[data.Rainfall .> 0], initialvalues = [log(5), 0.15, log(10)])

        # Parameter estimates
        θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]]

        # Parameter estimates with Naveau (R)
        θ = [5.30852235690306; 0.167910382987591; 29.9938206633414]

        @test θ̂ ≈ θ rtol = 0.01
    end

    @testset "Model (iv)" begin

        data = Extremes.dataset("rain")

        # Fit of the first model by maximum likelihood
        fm = EGPbetapowerfit(data.Rainfall[data.Rainfall .> 0], initialvalues = [log(5), 0.15, log(1.2), log(2)])

        # Parameter estimates
        θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]; params(fm)[4]]

        # Parameter estimates with Naveau (R)
        θ = [5.93923996559829; 0.12214112406368; 25.0590202053487; 1.80604091227634]

        @test θ̂ ≈ θ rtol = 0.01
    end
end
