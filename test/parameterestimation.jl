# @testset "parameterestimation" begin
#     @testset "Model (i)" begin

#         model = EGPpower(1, .3, 2)
#         data = rand(model, 10000)

#         # Fit of the first model by maximum likelihood
#         fm = EGPpowerfit(data)

#         # Parameter estimates
#         θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]]

#         # True parameters
#         θ = [params(model)[1]; params(model)[2]; params(model)[3]]

#         @test θ̂ ≈ θ rtol = 0.2
#     end

#     @testset "Model (iii)" begin

#         model = EGPbeta(1, .3, 2)
#         data = rand(model, 10000)

#         # Fit of the first model by maximum likelihood
#         fm = EGPbetafit(data)

#         # Parameter estimates
#         θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]]

#         # True parameters
#         θ = [params(model)[1]; params(model)[2]; params(model)[3]]

#         @test θ̂ ≈ θ rtol = 0.2
#     end

#     @testset "Model (iv)" begin

#         model = EGPbetapower(1, .3, 2,2)
#         data = rand(model, 10000)

#         # Fit of the first model by maximum likelihood
#         fm = EGPbetapowerfit(data, initialvalues = [5, 0.15, 2, 2])

#         # Parameter estimates
#         θ̂ = [params(fm)[1]; params(fm)[2]; params(fm)[3]; params(fm)[4]]

#         # True parameters
#         θ = [params(model)[1]; params(model)[2]; params(model)[3]; params(model)[4]]

#         @test θ̂ ≈ θ rtol = 0.2
#     end
# end
