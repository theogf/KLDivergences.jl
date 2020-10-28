# KLDivergences
[![Build Status](https://travis-ci.com/theogf/KLDivergences.jl.svg?branch=main)](https://travis-ci.com/theogf/KLDivergences.jl)
[![Coverage Status](https://coveralls.io/repos/github/theogf/KLDivergences.jl/badge.svg?branch=main)](https://coveralls.io/github/theogf/KLDivergences.jl?branch=main)


Compute [KL Divergences](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) using [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) objects either analytically (PR are welcome to enrich the library) or via Monte Carlo sampling.

There is only one function exported : `KL` which takes two arguments `p` and `q` for an analytical formulation and an additional argument `n_samples` for the sampling based approach.

The list of the pair of distributions supported analytically is given here 

## Univariate distributions

| p | q |
|:-------:|:-------:|
| Normal | Normal |
| Poisson | Poisson |
| Exponential | Exponential |
| Gamma | Gamma |
| InverseGamma | InverseGamma |
| Beta | Beta |

## Multivariate distributions

| p | q |
|:---:|:---:|
| MvNormal | MvNormal |
