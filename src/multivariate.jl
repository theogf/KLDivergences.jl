function KL(p::MvNormal, q::MvNormal)
    length(p) == length(q) || 
        throw(DimensionMismatch("Distributions p and q have different dimensions $(length(p)) and $(length(q))"))
    Σp = p.Σ; Σq = q.Σ
    0.5 * (tr(Σq \ Σp) + invquad(Σq, mean(p) - mean(q)) - length(p) + logdet(Σq) - logdet(Σp))
end