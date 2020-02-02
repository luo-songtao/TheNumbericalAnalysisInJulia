
"""
# 共轭梯度方法
    conjugate_gradient(A, b, x_0, n)

矩阵``A``是对称正定``n\\times n``矩阵，``b \\neq 0``是一个向量。

## 共轭梯度方法公式：
- 令``\\alpha_k = \\frac{r^T_kr_k}{d^T_kAd_k}``
- 令``\\beta_k = \\frac{r^T_{k+1}r_{k+1}}{r^T_kr_k}``
``\\begin{aligned} x_0 &= 初始估计 \\\\ d_0 &= r_0 = b-Ax_0 \\\\ x_{k+1} &= x_k + \\alpha_k d_k \\\\ r_{k+1} &= r_k - \\alpha_k A d_k \\\\ d_{k+1} &= r_{k+1} + \\beta_k d_k \\end{aligned}``
- 当``r_k=0``时，方程便得以求解

共轭梯度方法需要在每一步更新三个向量：
- ``x_k``: 表示第``k``步的近似解。
    - 更新公式：``x_{k+1} = x_k + \\alpha_k d_k```
- ``r_k``: 表示近似解``x_k``的余项。
    - 更新公式：``r_{k+1} = r_k - \\alpha_k A d_k``
- ``d_k``: 表示用于更新``x_k``得到改进的``x_{k+1}``时，所使用的新的搜索方向。
    - 更新公式：``d_{k+1} = r_{k+1} + \\beta_k d_k``

共轭梯度方法成功的关键：
- 所有的余项``r_k``两两正交: ``r^T_kr_j = 0,j<k``
- 方向``d_k``两两``A``共轭: ``d^T_k A d_j = 0,j<k``

同时为了保证下一个余项向量与前面所有余项向量都正交
- 余项``r_{k+1}``和方向``d_k``是正交的
- ``r_j``和方向``d_k``是A共轭的(j<k)

## 关于``\\alpha_k``的推导:
- 需要精确``\\alpha_k``使得新的余项``r_{k+1}``和方向``d_k``正交
``\\begin{aligned} x_{k+1} &= x_k + \\alpha_k d_k \\\\ b-Ax_{k+1} &= b-Ax_k - \\alpha_k A d_k \\\\ r_{k+1} &= r_k - \\alpha_k A d_k \\\\ 0 = d^T_kr_{k+1} &= d^T_k r_k- \\alpha_k d^T_k A d_k \\\\ \\alpha_k &= \\frac{d^T_kr_k}{d^T_kAd_k} \\end{aligned}``
- 因为``r_k``和方向``d_{k-1}``是正交
``\\begin{aligned} d_k-r_k &= \\beta_{k-1}d_{k-1} \\\\ r^T_kd_k - r^T_kr_k &= 0 \\end{aligned}``
- 可得``\\alpha_k = \\frac{r^T_kr_k}{d^T_kAd_k}``

## 关于``\\beta_k``的推导:
``\\begin{aligned} d_{k+1} &= r_{k+1} + \\beta_k d_k \\\\ 0 = d^T_k A d_{k+1} &= d^T_kAr_{k+1} +\\beta_k d^T_kAd_k \\\\ \\beta_k &= -\\frac{d^T_kAr_{k+1}}{d^T_kAd_k} \\end{aligned}``
- 因为余项向量两两正交
``\\begin{aligned} r_{k+1} &= r_k - \\alpha_k A d_k \\\\ r^T_jr_{k+1} &= r^T_jr_k - \\alpha_k r^T_j A d_k = 0 \\\\ \\frac{r^T_jr_k}{r^T_kr_k} &= - \\frac{r^T_j A d_k}{d^T_kAd_k} \\end{aligned}``
- 对于``j=k+1``
``\\beta_k = -\\frac{d^T_kAr_{k+1}}{d^T_kAd_k} = \\frac{r^T_{k+1}r_{k+1}}{r^T_kr_k}``

# Example
```jldoctest
julia> A = [2 2; 2 5]
2×2 Array{Int64,2}:
 2  2
 2  5
julia> b = [6;3]
2-element Array{Int64,1}:
 6
 3
julia> x = conjugate_gradient(A, b, [0;0], 2)
2-element Array{Float64,1}:
  4.0
 -0.9999999999999993
```
"""
function conjugate_gradient(A, b, x_0, n)
    x_k = x_0
    r_k = b - A*x_k    # 这里是r_0
    d_k = r_k    # 这里是d_0
    for k = 0:n-1
        if iszero(r_k); break; end 
        α_k = (transpose(r_k) * r_k)/(transpose(d_k)*A*d_k)
        x_k = x_k + α_k * d_k
        r_k_plus_1 = r_k - α_k* A * d_k
        β_k = (transpose(r_k_plus_1)*r_k_plus_1)/(transpose(r_k)*r_k)
        d_k = r_k_plus_1 + β_k * d_k
        r_k = r_k_plus_1
    end
    return x_k
end

"""
A = [1 -1 0; -1 2 1; 0 1 2]
b = [0;2;3]
x = conjugate_gradient(A, b, [0;0;0], 3)
println(x) # [0.9999999999999993, 1.0000000000000004, 1.0000000000000004]

A = [1 -1 0; -1 2 1; 0 1 5]
b = [3;-3;4]
x = conjugate_gradient(A, b, [0;0;0], 3)
println(x) # [1.9999999999999996, -0.9999999999999987, 1.0000000000000004]
"""