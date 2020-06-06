struct GewekeTest <: AbstractMCMCTest
    n_samples::Int
end

mutable struct GewekeTestResult{T} <: AbstractMCMCResult
    samples_fwd::T
    samples_bwd::T
    statistic
end

function Base.show(io::IO, res::GewekeTestResult)
    println(io, "Geweke (Joint Distribution) Test")
    println(io, "--------------------------------")
    println(io, "Results:")
    println(io, "    Number of samples: $(size(res.samples_fwd, 2))")
    println(io, "    Parameter dimension: $(size(res.samples_fwd.θ, 1))")
    println(io, "    Data dimension: $(size(res.samples_fwd.x, 1))")
    println(io, "    Statistic: $(res.statistic)")
    if ismissing(res.statistic)
        println(io, "")
        println(
            io, 
            """
            Test statistic is missing. Please use `compute_statistic!(res, g)` 
            if you want compute statistic without rerun the simulation.
            """
        )
    end
end

"""
    perform(cfg::GewekeTest, rand_θ, rand_x_given, rand_θ_given, g=nothing)

Run Geweke (joint distribution) test and compute the test statistic 
using `g` as the test function as in Equation (6) of (Geweke, 2014).
"""
function perform(cfg::GewekeTest, rand_θ, rand_x_given, rand_θ_given, g=nothing)
    @unpack n_samples = cfg
    
    # Generate samples
    local dim_θ, dim_x, samples_fwd, samples_bwd
    θ_bwd = rand_θ()
    @showprogress for i in 1:n_samples
        # Marginal-cnditional simulator
        θ_fwd = rand_θ()
        x_fwd = rand_x_given(θ_fwd)
        if i == 1
            dim_θ = length(θ_fwd)
            dim_x = length(x_fwd)
            dim = dim_θ + dim_x
            T = eltype(θ_fwd)
            samples_fwd = Matrix{T}(undef, dim, n_samples)
            samples_bwd = Matrix{T}(undef, dim, n_samples)
        end
        samples_fwd[:,i] = cat(θ_fwd, x_fwd; dims=1)
        # Successive-conditional simulator
        x_bwd = rand_x_given(θ_bwd)
        if i > 1
            θ_bwd = rand_θ_given(x_bwd)
        end
        samples_bwd[:,i] = cat(θ_bwd, x_bwd; dims=1)
    end
    samples_fwd = @LArray samples_fwd (θ=(1:dim_θ,:), x=(dim_θ+1:dim_θ+dim_x,:))
    samples_bwd = @LArray samples_bwd (θ=(1:dim_θ,:), x=(dim_θ+1:dim_θ+dim_x,:))

    # Compute statistics
    if isnothing(g)
        @warn "Test function `g` is not provided. Statistic is not computed."
        statistic = missing
    else
        statistic = _compute_statistic(samples_fwd, samples_bwd, g)
    end
    return GewekeTestResult(samples_fwd, samples_bwd, statistic)
end

function _compute_statistic(samples_fwd, samples_bwd, g)
    g_fwd = map(i -> g(samples_fwd.θ[:,i], samples_fwd.x[:,i]), 1:size(samples_fwd, 2))
    g_bwd = map(i -> g(samples_bwd.θ[:,i], samples_bwd.x[:,i]), 1:size(samples_bwd, 2))
    m_fwd = mean(g_fwd)
    v_fwd = mean(x -> x.^2, g_fwd) - m_fwd.^2
    m_bwd = mean(g_bwd)
    v_bwd = mean(x -> x.^2, g_bwd) - m_bwd.^2
    M₁, M₂ = length(samples_fwd), length(samples_bwd)
    statistic = (m_fwd - m_bwd) ./ sqrt.(v_fwd / M₁ + v_bwd / M₂)
end

function compute_statistic!(res::GewekeTestResult, g)
    @unpack samples_fwd, samples_bwd = res
    res.statistic = _compute_statistic(samples_fwd, samples_bwd, g)
    return res
end

mmd_of(res::MCMCDebugging.GewekeTestResult; kwargs...) = 
    mmd_of(res.samples_fwd, res.samples_bwd; kwargs...)

@recipe function f(res::GewekeTestResult, logjoint; n_grids=100)
    @unpack samples_fwd, samples_bwd = res
    logjoint_fwd = map(i -> logjoint(samples_fwd.θ[:,i], samples_fwd.x[:,i]), 1:size(samples_fwd, 2))
    logjoint_bwd = map(i -> logjoint(samples_bwd.θ[:,i], samples_bwd.x[:,i]), 1:size(samples_bwd, 2))
    joint_fwd, joint_bwd = exp.(logjoint_fwd), exp.(logjoint_bwd)

    # Compute CDFs
    mins, maxs = extrema.([joint_fwd, joint_bwd])
    percent_fwd = Vector{Float64}(undef, n_grids)
    percent_bwd = Vector{Float64}(undef, n_grids)
    for (i, v) in enumerate(range(minimum(mins), maximum(maxs); length=n_grids))
        percent_fwd[i] = sum(joint_fwd .< v) / length(joint_fwd)
        percent_bwd[i] = sum(joint_bwd .< v) / length(joint_bwd)
    end
    
    # Recipe
    legend --> :topleft
    
    @series begin
        linecolor := 2
        label := "Sampler"
        percent_fwd, percent_bwd
    end
    
    @series begin
        linecolor := :gray
        linestyle := :dash
        label := "Perfect"
        [0, 1], [0, 1]
    end
end