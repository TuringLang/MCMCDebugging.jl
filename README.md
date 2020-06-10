# MCMCDebugging.jl: debugging utilities for MCMC samplers

This package implements a few utilities for debugging MCMC samplers, which includes

- Geweke test
  - See the references [1,2] or [this blog](https://lips.cs.princeton.edu/testing-mcmc-code-part-2-integration-tests/) for details
- [ ] Central limit theorem test

See the [notebook](docs/example.ipynb) for an example.

## TODOs

- [ ] Interface with [DynamicPPL.jl](https://github.com/TuringLang/DynamicPPL.jl) so that `rand_θ` and `rand_x_given` can be automatically generated.
- [ ] Interface with [AbstractMCMC.jl](https://github.com/TuringLang/AbstractMCMC.jl) so that `rand_θ_given` can be automatically generated.
- [ ] Support KSD for Geweke test via [KernelGoodnessOfFit.jl](https://github.com/torfjelde/KernelGoodnessOfFit.jl/tree/master/src).
  - MMD seems to work quite well so I would expect KSD to work better.

## References

[1] Geweke J. Getting it right: Joint distribution tests of posterior simulators. Journal of the American Statistical Association. 2004 Sep 1;99(467):799-804.

[2] Grosse RB, Duvenaud DK. Testing mcmc code. arXiv preprint arXiv:1412.5218. 2014 Dec 16.
