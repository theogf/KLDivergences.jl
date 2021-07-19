module KLDivergences

using Distributions: StatsBase
using Distributions
using LinearAlgebra
using PDMats
using Distances
using SpecialFunctions
using StatsBase: StatsBase, kldivergence


export KL, kldivergence

"""
    KL(p::Distribution, q::Distribution) -> T
    KL(p::Distribution, q::Distribution, n_samples::Int) -> T

Return the KL divergence of  KL(p||q), either by sampling or analytically
"""
KL

StatsBase.kldivergence(p::Sampleable, q::Sampleable) = KL(p, q)

KLbase(p, q, x) = logpdf(p, x) - logpdf(q, x)

## Generic fallback for multivariate Distributions
function KL(p::UnivariateDistribution, q::UnivariateDistribution, n_samples = 1_000)
    return mean(x->KLbase(p, q, x), rand(p, n_samples))
end

function KL(p::MultivariateDistribution, q::MultivariateDistribution, n_samples = 1_000)
    length(p) == length(q) ||
        throw(DimensionMismatch("Dimensions of p and q do not match"))
    return mean(x->KLbase(p, q, x), eachcol(rand(p, n_samples)))
end

function symmetricKL(p, q)
    return KL(p, q) + KL(q, p)
end

include("univariate.jl")
include("multivariate.jl")

end
