using Distributions, DynamicPPL

# The marginal-conditional simulator defined by DynamicPPL
# See README.md for the expected form of the model definition.
@model function BetaBinomial(θ=missing, x=missing)
    θ ~ Beta(2, 3)
    x ~ Binomial(3, θ)
    return θ, x
end

# The successive-conditional simulator
# 1. Bug-free posterior sampler
# Beta(α + x, β + n - x) is the true posterior.
rand_θ_given(x) = rand(Beta(2 + x, 3 + 3 - x))
# 2. Buggy posterior sampler
rand_θ_given_buggy(x) = rand_θ_given(min(3, x + 1))

# Test function
g(θ, x) = cat(θ, x; dims=1)

using MCMCDebugging

res = perform(GewekeTest(5_000), BetaBinomial, rand_θ_given; g=g)

using Plots

plot(res, BetaBinomial(); size=(300, 300), title="Bug-free sampler")

res_buggy = perform(GewekeTest(5_000), BetaBinomial, rand_θ_given_buggy)

compute_statistic!(res_buggy, g)

plot(res_buggy, BetaBinomial(); size=(300, 300), title="Buggy sampler")

@info "MMD" mmd_of(res) mmd_of(res_buggy)
