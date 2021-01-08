module MCMCDebugging

using UnPack, ProgressMeter, RecipesBase, Statistics, LabelledArrays, 
    Distributions, HypothesisTests, DynamicPPL

abstract type AbstractMCMCTest end
abstract type AbstractMCMCResult end

include("mmd.jl")
export mmd_of

include("geweke.jl")
export GewekeTest

include("clt.jl")
export CLTTest

### DynamicPPL integration

function perform(cfg::GewekeTest, modelgen::Function, rand_θ_given::Function; kwargs...)
    model = modelgen()
    rand_marginal() = model()
    function rand_x_given(θ)
        _, x = modelgen(θ)()
        return x
    end
    return perform(cfg, rand_marginal, rand_x_given, rand_θ_given; kwargs...)
end

@recipe function f(res::GewekeTestResult, model::Model)
    vi = VarInfo(model)
    spl = SampleFromPrior()
    function _logjoint(θ, x)
        vi[spl] = cat(θ, x; dims=1)
        return logjoint(model, vi)
    end
    res, _logjoint
end

export perform, compute_statistic!

end # module
