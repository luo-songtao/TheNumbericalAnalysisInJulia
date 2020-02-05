"""
# Broyden方法
    broyden(F, A_0, x_0, k)
当无法得到雅可比矩阵``DF``时，Broyden方法通常被认为是一种次优的方法。

将设A_i是第i步可以得到的雅可比矩阵的最优近似，并被用于：
```math
x_{k+1} = x_k - A^{-1}_iF(x_i)
```

为了根据``A_i``更新得到``A_{i+1}``。首先注意到雅可比矩阵满足：

``\\begin{aligned} &A_{i+1}\\delta_{i+1} = \\Delta_{i+1} \\\\ &\\delta_{i+1} = x_{i+1} - x_i \\\\ &\\Delta_{i+1} = F(x_{i+1}) - F(x_i) \\end{aligned}``

要求所有满足``\\delta^T_{i+1}\\omega = 0``的``\\omega``：

```math
A_{i+1}\\omega = A_i \\omega
```

则同时满足上述条件的矩阵如下:

```math
A_{i+1} = A_i + \\frac{(\\Delta_{i+1} - A_i\\delta_i)\\delta^T_{i+1}}{\\delta^T_{i+1}\\delta_{i+1}}
```

# Example
```jldoctest
# F = [(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)]
# A_0 = Matrix{Float64}(I, 2, 2)
# x_0 = Float64[1, 2]
julia> using LinearAlgebra

julia> broyden([(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)],  Matrix{Float64}(I, 2, 2), Float64[1, 2], 20)
2×1 Array{Float64,2}:
 0.8260271593976959
 0.5635461739578368
```
"""
function broyden(F, A_0, x_0, k)
    n = length(x_0)
    x_k = x_0
    A_k = A_0
    y_k = Float64[f(x_k...) for f in F]
    for i = 1:k
        s = naive_gauss_elimination(A_k, y_k, n)
        x_k = x_k - s
        δ_i_plus_1 = - s
        y_k_plus_1 = Float64[f(x_k...) for f in F]
        Δ_i_plus_1 = y_k_plus_1 - y_k
        A_k = A_k + (Δ_i_plus_1 - A_k*δ_i_plus_1)*transpose(δ_i_plus_1) ./ (transpose(δ_i_plus_1) * δ_i_plus_1)
        y_k = y_k_plus_1
    end
    return x_k
end

include("naive_gauss_elimination.jl")
