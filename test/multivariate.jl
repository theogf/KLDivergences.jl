@testset "univariate" begin
    struct CholeskyMvNormal{TL,Tm} <: Distributions.AbstractMvNormal
        m::Tm
        L::TL
    end
    Distributions.mean(p::CholeskyMvNormal) = p.m
    Distributions.cov(p::CholeskyMvNormal) = p.L * p.L'
    Distributions.rand(p::CholeskyMvNormal, n::Int) = p.m .+ p.L * randn(length(p), n) 
    Distributions.length(p::CholeskyMvNormal) = length(p.m)
    function Distributions.logpdf(p::CholeskyMvNormal, x::AbstractVector)
        return -0.5 * (length(p) * log(2π) + 2 * logdet(p.L) + sum(abs2, p.L \ (x .- p.m)))
    end
    @testset "AbstractMvNormal" begin
        n_dim = 2
        X1 = cholesky(Matrix(0.5 * I(n_dim))).L
        X2 = cholesky(Matrix(0.3 * I(n_dim))).L
        p = CholeskyMvNormal(zeros(n_dim), X1)
        q = CholeskyMvNormal(ones(n_dim), X2)
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
    @testset "MvNormal" begin
        n_dim = 2
        p = MvNormal(zeros(n_dim), Matrix(0.5 * I(n_dim)))
        q = MvNormal(ones(n_dim), Matrix(0.3 * I(n_dim)))
        @test KL(p, q) > 0
        @test KL(p, q) ≈ KL(p, q, 100_000) atol = 0.1
    end
end