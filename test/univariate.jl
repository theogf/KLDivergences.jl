@testset "univariate" begin
    @testset "Beta" begin
        p = Beta(2, 10)
        q = Beta(3, 5)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "Exponential" begin
        p = Exponential(2.0)
        q = Exponential(3.0)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "Gamma" begin
        p = Gamma(2.0, 1.0)
        q = Gamma(3.0, 2.0)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "InverseGamma" begin
        p = InverseGamma(2.0, 1.0)
        q = InverseGamma(3.0, 2.0)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "Normal" begin
        p = Normal(0, 1)
        q = Normal(0.5, 0.5)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "Poisson" begin
        p = Poisson(4.0)
        q = Poisson(3.0)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    
end