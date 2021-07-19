using KLDivergences
using Distributions
using LinearAlgebra
using Random
using Test

@testset "KLDivergences.jl" begin
    include("univariate.jl")
    include("multivariate.jl")
end