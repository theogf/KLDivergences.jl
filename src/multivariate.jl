function KL(p::AbstractMvNormal, q::AbstractMvNormal)
    length(p) == length(q) || 
        throw(DimensionMismatch("Distributions p and q have different dimensions $(length(p)) and $(length(q))"))
    Σp = cov(p)
    Σq = cov(q)
    0.5 * (tr(Σq \ Σp) + invquad(Σq, mean(p) - mean(q)) - length(p) + logdet(Σq) - logdet(Σp))
end