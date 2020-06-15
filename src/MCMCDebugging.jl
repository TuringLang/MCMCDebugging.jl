module MCMCDebugging

using Parameters, ProgressMeter, RecipesBase, Statistics, LabelledArrays, 
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

perform(cfg::GewekeTest, modelgen::ModelGen, rand_θ_given, g=nothing) = 
    perform(cfg, modelgen(missing, missing), θ -> modelgen(θ, missing)()[2], rand_θ_given, g)

@recipe function f(res::GewekeTestResult, modelgen::ModelGen)
    m = modelgen(missing, missing)
    vi = VarInfo(m)
    spl = SampleFromPrior()
    function _logjoint(θ, x)
        vi[spl] = cat(θ, x; dims=1)
        return logjoint(m, vi)
    end
    res, _logjoint
end

export perform, compute_statistic!

end # module
