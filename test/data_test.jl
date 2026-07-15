@testset "dataset()" begin
    import ExtendedExtremes.dataset

    for name in ("pcp", "tasmax")
        data = dataset(name)

        @test data isa DataFrame
        @test !isempty(data)
    end

    @test_throws ArgumentError dataset("unknown")

end