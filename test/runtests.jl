using Distributions: kldivergence
using KLDivergences
using Distributions
using LinearAlgebra
using Random
using StatsBase
using Test

@testset "KLDivergences.jl" begin
    include("univariate.jl")
    include("multivariate.jl")

    @testset "Generic Methods" begin
        p = Exponential(2.0)
        q = Exponential(5.0)
        @test symmetricKL(p, q) == symmetricKL(q, p)
        @test kldivergence(p, q) == KL(p, q)
    end

end