# This file is currently required to get AD on KL(::MvNormal, ::MvNormal) to work

import ChainRulesCore

# opt-out of generic rrule for `Matrix` defined in ChainRules.jl.
# It's horrible that this code is needed. Given that it does, its ideal home is PDMats.jl.
# We're type-pirating for now.
ChainRulesCore.@opt_out ChainRulesCore.rrule(::Type{<:Matrix}, ::PDMats.PDMat)
