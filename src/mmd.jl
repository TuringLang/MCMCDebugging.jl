function euclidsq(X::T, Y::T) where {T<:AbstractMatrix}
    XiXj = transpose(X) * Y
    x² = sum(X .^ 2; dims=1)
    y² = sum(Y .^ 2; dims=1)
    transpose(x²) .+ y² - 2XiXj
end

function euclidsq(X::T) where {T<:AbstractMatrix}
    XiXj = transpose(X) * X
    x² = sum(X .^ 2; dims=1)
    transpose(x²) .+ x² - 2XiXj
end

gaussian_gramian(esq, σ::AbstractFloat) = exp.(-esq ./ 2σ^2)

"""
    mmd_of(x_nu, x_de; σ=nothing)

Compute the maximum mean discrepency give two set of samples `x_nu` and `x_de`.
If `σ` is not provided, it will be taken as "median / log(n)" where `median` is 
the median of pair-wise Ecludean distances and `n` is the number of samples in `x_de`.
"""
function mmd_of(x_nu, x_de; σ=nothing)
    d²_dede, d²_denu, d²_nunu = euclidsq(x_de), euclidsq(x_de, x_nu), euclidsq(x_nu)
    # Heuristic: take `σ²` as "median / log(n)"
    if σ isa Nothing
        h = median(vcat(vec.([d²_dede, d²_denu, d²_nunu]))) / log(size(d²_dede, 1))
        σ = sqrt(h)
    end
    Kdede, Kdenu, Knunu = gaussian_gramian.((d²_dede, d²_denu, d²_nunu), σ)
    return mean(Kdede) - 2mean(Kdenu) + mean(Knunu)
end
