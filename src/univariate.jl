"""
    KL(p::Beta, q::Beta)

See [KL Beta](https://en.wikipedia.org/wiki/Beta_distribution#Quantities_of_information_(entropy))
"""
function KL(p::Beta, q::Beta)
    αp, βp = params(p)
    αq, βq = params(q)
    return logbeta(αq, βq) - logbeta(αp, βp) + (αp - αq) * digamma(αp) +
        (βp - βq) * digamma(βp) + (αq - αp + βq - βp) * digamma(αp + βp)
end

"""
    KL(p::Exponential, q::Exponential)

See [KL Exponential](https://en.wikipedia.org/wiki/Exponential_distribution#Kullback%E2%80%93Leibler_divergence)
"""
function KL(p::Exponential, q::Exponential)
    λp = scale(p)
    λq = scale(q)
    return log(λp) - log(λq) + λq / λp - 1
end

"""
    KL(p::Gamma, q::Gamma)

See [KL Gamma](https://en.wikipedia.org/wiki/Gamma_distribution#Kullback%E2%80%93Leibler_divergence)
"""
function KL(p::Gamma, q::Gamma)
    # We use the parametrization with the rate β
    αp, αq = shape.((p, q))
    βp, βq = rate.((p, q))
    return (αp - αq) * digamma(αp) - loggamma(αp) + loggamma(αq) +
        αq * (log(βp) - log(βq)) + αp * (βq - βp) / βp
end

"""
    KL(p::InverseGamma, q::InverseGamma)

See [KL Inverse-Gamma](https://en.wikipedia.org/wiki/Inverse-gamma_distribution#Properties)
"""
function KL(p::InverseGamma, q::InverseGamma)
    # We can reuse the implementation of Gamma
    return KL(Gamma(shape(p), rate(p)), Gamma(shape(q), rate(q)))
end

"""
    KL(p::Normal, q::Normal)

See [KL Gaussian](https://en.wikipedia.org/wiki/Normal_distribution#Other_properties)
"""
function KL(p::Normal, q::Normal)
    μp, σp = params(p)
    μq, σq = params(q)
    return 0.5 * (abs2(σp / σq) + abs2((μp - μq) / σq) - 1 + 2 * (log(σq) - log(σp)))
end

"""
    KL(p::Poisson, q::Poisson)

See [KL Poisson](https://en.wikipedia.org/wiki/Poisson_distribution#Other_properties)
"""
function KL(p::Poisson, q::Poisson)
    λp, λq = rate.((p, q))
    return λq - λp + λp * (log(λp) - log(λq))
end



