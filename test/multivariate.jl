@testset "univariate" begin
    @testset "MvNormal" begin
        n_dim = 2
        p = MvNormal(zeros(n_dim), Matrix(0.5 * I(n_dim)))
        q = MvNormal(ones(n_dim), Matrix(0.3 * I(n_dim)))
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 1.0
    end
end