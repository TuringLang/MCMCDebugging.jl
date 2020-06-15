# MCMCDebugging.jl: debugging utilities for MCMC samplers

This package implements a few utilities for debugging MCMC samplers, which includes

- [x] Geweke test
  - See the references [1,2] or [this blog](https://lips.cs.princeton.edu/testing-mcmc-code-part-2-integration-tests/) for details
- [ ] Central limit theorem test

See the [notebook](https://nbviewer.jupyter.org/github/xukai92/MCMCDebugging.jl/blob/master/docs/example.ipynb) for an example.

## Usage

### Defining test models via DynamicPPL.jl

MCMCDebugging.jl allows using DynamicPPL.jl to define test models.
In the example [notebook](https://nbviewer.jupyter.org/github/xukai92/MCMCDebugging.jl/blob/master/docs/example.ipynb), the test model is defined as

```julia
@model function BetaBinomial(θ, x)
    θ ~ Beta(2, 3)
    x ~ Binomial(3, θ)
    return θ, x
end
```

There are a few requirements from MCMCDebugging.jl to use the defined model.

1. The model should take `θ` and `x` as inputs.
2. The model should return the parameter `θ` and the data `x` as a tuple.

With these two points, MCMCDebugging.jl can generate several functions used by lower-level APIs.

1. `rand_marginal()`: drawing `θ` and `x` as a tuple
2. `rand_x_given(θ)`: drawing `x` conditioned on `θ`
3. `logjoint(θ, x)`: computing the log-joint probability of `θ` and `x`

1 and 2 are used to perform the Geweke test and 3 is used to make the Q-Q plot.

## Lower-level APIs

### Geweke test

Defining the Geweke test

```julia
cfg = GewekeTest(n_samples::Int)
```

where `n_samples` is the number of samples used for testing.

Performing the Geweke test

```julia
res = perform(cfg::GewekeTest, rand_marginal, rand_x_given, rand_θ_given, g=nothing)
```

where

- `rand_marginal()` draws `θ` and `x` as a tuple
- `rand_x_given(θ)` draws `x` conditioned on `θ`
- `rand_θ_given(x)` draws `θ` conditioned on `x`
- `g(θ, x)` is the test function

Making the Q-Q plot

```julia
plot(res::GewekeTestResult, logjoint)
```

where

- `logjoint(θ, x)` computes the log-joint probability of `θ` and `x`

## TODOs

- [x] Interface with [DynamicPPL.jl](https://github.com/TuringLang/DynamicPPL.jl) so that `rand_marginal` and `rand_x_given` can be automatically generated.
- [ ] Interface with [AbstractMCMC.jl](https://github.com/TuringLang/AbstractMCMC.jl) so that `rand_θ_given` can be automatically generated.
- [ ] Support KSD for Geweke test via [KernelGoodnessOfFit.jl](https://github.com/torfjelde/KernelGoodnessOfFit.jl/tree/master/src).

## References

[1] Geweke J. Getting it right: Joint distribution tests of posterior simulators. Journal of the American Statistical Association. 2004 Sep 1;99(467):799-804.

[2] Grosse RB, Duvenaud DK. Testing mcmc code. arXiv preprint arXiv:1412.5218. 2014 Dec 16.
