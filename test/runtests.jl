using KLDivergences
using Distributions
using LinearAlgebra
using Test

@testset "KLDivergences.jl" begin
    include("univariate.jl")
    include("multivariate.jl")
    # Write your tests here.
end