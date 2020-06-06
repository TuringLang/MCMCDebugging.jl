struct CLTTest <: AbstractMCMCTest
    n_chains::Int
end

struct CLTTestResult{T} <: AbstractMCMCResult
    chains::T
end
