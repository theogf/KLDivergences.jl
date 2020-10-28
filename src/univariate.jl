# λq - λp + λp log λp / λq
function KL(p::Poisson, q::Poisson)
    mean(q) - mean(p) + mean(p) * (log(mean(p)) - log(mean(q)))
end

function KL(p::Normal, q::Normal)
    0.5 * (var(p) / var(q) + abs2(mean(p) - mean(q)) / var(q) - 1 + 2 * (log(std(q)) - log(std(p))))
end

function KL()

end