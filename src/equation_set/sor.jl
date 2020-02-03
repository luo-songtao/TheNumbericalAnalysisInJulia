"""
# 连续过松弛(SOR)
    sor(A, b, x_0, n, k, ω=1.25)

连续过松弛方法使用高斯-赛德尔方法的求解方向，并使用过松弛以加快收敛速度。

令``\\omega``是一个实数，将新的估计中的每个元素``x_{k+1}``定义为``\\omega``乘上高斯-赛德尔公式和``1-\\omega`乘上当前估计的平均。

``\\omega``被称为松弛参数，当``\\omega>1``时被称为过松弛。

当``\\omega=1``时，SOR方法就是高斯-赛德尔方法

## 连续过松弛(SOR)迭代公式：
``\\begin{aligned} x_0 &= 初始估计 \\\\ x_{k+1} &= (1-\\omega)x_k + \\omega D^{-1}(b-Lx_{k+1}-Ux_k), k=0,1,2,...... \\end{aligned}``

# Arguments
- `A`: 表示系数矩阵A
- `b`: 表示常数项b
- `x_0`: 初始估计(向量)
- `n`: 方程数
- `k`: 迭代次数
- `ω`: 松弛参数(default=1.25)

# Example
```jldoctest
julia> A = Float64[3 1 -1; 2 4 1; -1 2 5]
3×3 Array{Float64,2}:
  3.0  1.0  -1.0
  2.0  4.0   1.0
 -1.0  2.0   5.0
julia> b = Float64[4; 1; 1]
3-element Array{Float64,1}:
 4.0
 1.0
 1.0
julia> x_0 = Float64[0;0;0]
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
julia> x = sor(A, b, x_0, 3, 50, 1.25)    # 该例中约30步收敛
3-element Array{Float64,1}:
  2.0
 -1.0
  1.0
```
"""
function sor(A, b, x_0, n, k, ω=1.25)
    x_k = x_0
    for i = 1:k
        for j = 1:n
            L_mult = reshape(A[j, 1:j-1], (1, j-1)) * reshape(x_k[1:j-1], (j-1, 1))
            U_mult = reshape(A[j, j+1:n], (1, n-j)) * reshape(x_k[j+1:n], (n-j, 1))
            x_k[j] = (1-ω)x_k[j] + ω*((1/A[j,j]) * (b[j] - L_mult[1] - U_mult[1]))
        end
    end
    return x_k
end
