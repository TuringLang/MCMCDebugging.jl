# MCMCDebugging.jl: debugging utilities for MCMC samplers

This package implements a few utilities for debugging MCMC samplers, which includes

- [x] Geweke test
  - See the references [1,2] or [this blog](https://lips.cs.princeton.edu/testing-mcmc-code-part-2-integration-tests/) for details
- [ ] Central limit theorem test

See the [notebook](https://nbviewer.jupyter.org/github/xukai92/MCMCDebugging.jl/blob/master/docs/example.ipynb) for an example.

## Usage

The example [notebook](https://nbviewer.jupyter.org/github/xukai92/MCMCDebugging.jl/blob/master/docs/example.ipynb) covers most of the usages.
Some details on the model definition via DynamicPPL is explained below.

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
res = perform(cfg::GewekeTest, rand_marginal, rand_x_given, rand_θ_given, g=nothing; progress=true)
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

## References

[1] Geweke J. Getting it right: Joint distribution tests of posterior simulators. Journal of the American Statistical Association. 2004 Sep 1;99(467):799-804.

[2] Grosse RB, Duvenaud DK. Testing mcmc code. arXiv preprint arXiv:1412.5218. 2014 Dec 16.
