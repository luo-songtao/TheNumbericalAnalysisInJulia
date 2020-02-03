"""
# BroydenⅡ 方法

Broyden的第二种方式避免相对代价较大的矩阵求解步骤。相比近似``DF``导数，同样我们可以近似``DF``的逆

基于``B_i = A^{-1}_i``，我们希望得到：
```math
\\delta_{i+1} = B_{i+1}\\Delta_{i+1}
```

同样为了根据B_i更新得到B_{i+1}。首先注意到雅可比矩阵满足：
``\\begin{aligned} &A_{i+1}\\delta_{i+1} = \\Delta_{i+1} \\\\ &\\delta_{i+1} = x_{i+1} - x_i \\\\ &\\Delta_{i+1} = F(x_{i+1}) - F(x_i) \\end{aligned}``

要求所有满足``\\delta^T_{i+1}\\omega = 0``的``\\omega``：

```math
A_{i+1}\\omega = A_i \\omega
```
或者：

```math
B_{i+1}A_i\\omega = \\omega
```

则同时满足上述条件的矩阵如下:

```math
B_{i+1} = B_i + \\frac{(\\delta_{i+1} - B_i\\Delta_i)\\delta^T_{i+1}B_i}{\\delta^T_{i+1}B_i\\Delta_{i+1}}
```

不再需要进行矩阵求解的迭代公式：
```math
x_{i+1} = x_i - B_iF(x_i)
```

# Example
```jldoctest
# F = [(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)]
# B_0 = Matrix{Float64}(I, 2, 2)
# x_0 = Float64[1, 2]
julia> using LinearAlgebra

julia> broyden2([(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)],  Matrix{Float64}(I, 2, 2), Float64[1, 2], 50)
2-element Array{Float64,1}:
 0.8260313576541872
 0.5636241621612585
```
"""
function broyden2(F, B_0, x_0, k)
    x_k = x_0
    B_k = B_0
    y_k = Float64[f(x_k...) for f in F]
    for i = 1:k
        x_k_plus_1 = x_k - B_k * y_k
        δ_i_plus_1 = x_k_plus_1 - x_k
        y_k_plus_1 = Float64[f(x_k_plus_1...) for f in F]
        Δ_i_plus_1 = y_k_plus_1 - y_k

        B_k = B_k + (δ_i_plus_1 - B_k*Δ_i_plus_1)*transpose(δ_i_plus_1)*B_k ./ (transpose(δ_i_plus_1) * B_k * Δ_i_plus_1)
        y_k = y_k_plus_1
        x_k = x_k_plus_1

    end
    return x_k
end
