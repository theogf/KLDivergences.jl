@testset "univariate" begin
    @testset "MvNormal" begin
        n_dim = 2
        p = MvNormal(zeros(n_dim), rand(n_dim))
        q = MvNormal(ones(n_dim), rand(n_dim))
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 10_000) atol = 0.1
    end
end