function KL(p::Beta, q::Beta)
    return logbeta(q.α, q.β) - logbeta(q.α, q.β) + (p.α - q.α) * digamma(p.α) +
        (p.β - q.β) * digamma(p.β) + (q.α - p.α + q.β - p.β) * digamma(p.α + p.β)
end

function KL(p::Exponential, q::Exponential)
    return log(p.θ) - log(q.θ) + q.θ / p.θ - 1
end

function KL(p::Gamma, q::Gamma)
    return (p.α - q.α) * digamma(p.α) - loggamma(p.α) + loggamma(q.α) +
        q.α * (log(q.θ) - log(p.θ)) + p.α * (p.θ - q.θ) / q.θ
end

function KL(p::InverseGamma, q::InverseGamma)
    αp = p.invd.α; αq = q.invd.α
    return (αp - αq) * digamma(αp) - loggamma(αp) + loggamma(αq) +
        αq * (log(p.θ) - log(q.θ)) + αp * (q.θ - p.θ) / p.θ
end

function KL(p::Normal, q::Normal)
    return 0.5 * (var(p) / var(q) + abs2(mean(p) - mean(q)) / var(q) - 1 + 2 * (log(std(q)) - log(std(p))))
end

# λq - λp + λp log λp / λq
function KL(p::Poisson, q::Poisson)
    return q.λ - p.λ + p.λ * (log(p.λ) - log(q.λ))
end



