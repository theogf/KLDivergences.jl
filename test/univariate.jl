@testset "univariate" begin
    @testset "Normal" begin
        p = Normal(0, 1)
        q = Normal(0.5, 0.5)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 10_000) atol = 0.1
    end
    @testset "Poisson" begin
        p = Poisson(4)
        q = Normal(7.0)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 10_000) atol = 0.2
    end

end