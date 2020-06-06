module MCMCDebugging

using Parameters, ProgressMeter, RecipesBase, Statistics, LabelledArrays

abstract type AbstractMCMCTest end
abstract type AbstractMCMCResult end

include("mmd.jl")
export mmd_of

include("geweke.jl")
export GewekeTest

include("clt.jl")
export CLTTest

export perform, compute_statistic!

end # module
