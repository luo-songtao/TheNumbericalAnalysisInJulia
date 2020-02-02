"""
# 预条件共轭梯度方法
    conjugate_gradient_with_pre_condition(A, b, x_0, n)

共轭梯度法在病态矩阵上性能很差，可以通过**预条件**来得到缓解。主要是将问题转化为良态矩阵系统然后再实施共轭梯度法。

矩阵``A``是对称正定``n\\times n``矩阵，``b \\neq 0``是一个向量。

预条件形式：
```math
M^{-1}Ax = M^{-1}b
```
其中``M``是可逆的``n \\times n``矩阵，称为**预条件因子**。矩阵``M``应满足：
- 与矩阵``A``足够接近
- 容易求逆

但是和矩阵``A``最接近的是``A``自身，使用``M=A``会把问题条件数变为1，但是一般``A``不容易求逆。而最容易求逆的矩阵是单位矩阵，但它又不能减低条件数。

因此最完美的预条件矩阵位于两者的中间，同时具备二者的性质。
- 一种特别简单的方式是**雅可比预条件子``M=D``**（``D``是``A``的对角线矩阵）

## 预条件共轭梯度方法公式：
- 令``\\alpha_k = \\frac{r^T_kz_k}{d^T_kAd_k}``
- 令``\\beta_k = \\frac{r^T_{k+1}z_{k+1}}{r^T_kz_k}``
``\\begin{aligned} x_0 &= 初始估计 \\\\ d_0 &= r_0 = b-Ax_0 \\\\ x_{k+1} &= x_k + \\alpha_k d_k \\\\ r_{k+1} &= r_k - \\alpha_k A d_k \\\\ z_{k+1} &= M^{-1}r_{k+1} \\\\ d_{k+1} &= z_{k+1} + \\beta_k d_k \\end{aligned}``
- 当``r_k=0``时，``x_k``就是方程组的解

预条件共轭梯度方法需要在每一步更新四个向量：
- ``x_k``: 表示第``k``步的近似解。
    - 更新公式：``x_{k+1} = x_k + \\alpha_k d_k```
- ``r_k``: 表示近似解``x_k``的余项。
    - 更新公式：``r_{k+1} = r_k - \\alpha_k A d_k``
- ``z_k``: 表示预条件系统的余项。``z_k = M^{-1}b - M^{-1}Ax_k = m^{-1}r_k``
    - 更新公式：``z_{k+1} = M^{-1}r_{k+1}``
- ``d_k``: 表示用于更新``x_k``得到改进的``x_{k+1}``时，所使用的新的搜索方向。
    - 更新公式：``d_{k+1} = z_{k+1} + \\beta_k d_k``

# Example
```jldoctest
julia> using LinearAlgebra

julia> A = [1 -1 0; -1 2 1; 0 1 2]
3×3 Array{Int64,2}:
  1  -1  0
 -1   2  1
  0   1  2
julia> b = [0;2;3]
3-element Array{Int64,1}:
 0
 2
 3
julia> M = diagm(diag(A))    # 雅可比预条件子
3×3 Array{Int64,2}:
 1  0  0
 0  2  0
 0  0  2
julia> x = conjugate_gradient_with_pre_condition(A, b, [0;0;0], 3, M)
3-element Array{Float64,1}:
 0.9999999999999996
 0.9999999999999998
 0.9999999999999999

julia> A = [1 -1 0; -1 2 1; 0 1 5]
3×3 Array{Int64,2}:
  1  -1  0
 -1   2  1
  0   1  5
julia> b = [3;-3;4]
3-element Array{Int64,1}:
  3
 -3
  4
julia> M = diagm(diag(A))    # 雅可比预条件子
3×3 Array{Int64,2}:
 1  0  0
 0  2  0
 0  0  5
julia> x = conjugate_gradient_with_pre_condition(A, b, [0;0;0], 3, M)
3-element Array{Float64,1}:
  1.9999999999999996
 -1.0
  0.9999999999999999
```
"""
function conjugate_gradient_with_pre_condition(A, b, x_0, n, M)
    M_inverse = inv(M)
    x_k = x_0
    r_k = b - A*x_k    # 这里是r_0
    z_k = M_inverse*r_k    # 这里是z_0
    d_k = z_k    # 这里是d_0
    for k = 0:n-1
        if iszero(r_k); break; end 
        α_k = (transpose(r_k) * z_k)/(transpose(d_k)*A*d_k)
        x_k = x_k + α_k * d_k
        r_k_plus_1 = r_k - α_k* A * d_k
        z_k_plus_1 = M_inverse * r_k_plus_1
        β_k = (transpose(r_k_plus_1)*z_k_plus_1)/(transpose(r_k)*z_k)
        d_k = z_k_plus_1 + β_k * d_k
        r_k = r_k_plus_1
        z_k = z_k_plus_1
    end
    return x_k
end
