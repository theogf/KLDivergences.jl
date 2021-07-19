function KL(p::AbstractMvNormal, q::AbstractMvNormal)
    length(p) == length(q) || 
        throw(DimensionMismatch("Distributions p and q have different dimensions $(length(p)) and $(length(q))"))
    Σp = cov(p)
    Σq = cov(q)
    Δμ = mean(p) - mean(q)
    0.5 * (tr(Σq \ Σp) + dot(Δμ / Σq, Δμ) - length(p) + logdet(Σq) - logdet(Σp))
end

function KL(p::MvNormal, q::MvNormal)
    length(p) == length(q) || 
        throw(DimensionMismatch("Distributions p and q have different dimensions $(length(p)) and $(length(q))"))
    0.5 * (tr(q.Σ \ p.Σ) + invquad(q.Σ, mean(p) - mean(q)) - length(p) + logdet(q.Σ) - logdet(p.Σ))
end