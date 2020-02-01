
"""
# 楚列斯基分解
    cholesky_decomposition(A, n)

**楚列斯基分解定理**：``A``是``n \\times n``对称正定矩阵，则存在上三角矩阵``R``满足``A=R^TR``

``A = \\left [ \\begin{matrix} a & b \\\\ b & c \\end{matrix} \\right ] = \\left [ \\begin{matrix} \\sqrt{a} & 0 \\\\ \\frac{b}{\\sqrt{a}} & \\sqrt{c-b^2/a} \\end{matrix} \\right ] \\left [ \\begin{matrix} \\sqrt{a} & \\frac{b}{\\sqrt{a}} \\\\ 0 & \\sqrt{c-b^2/a} \\end{matrix} \\right ]=R^TR``

对于对称正定矩阵``A``，求解``Ax=b``，和``LU``分解方式相同：

``R^Tc=b \\Longrightarrow c \\\\ Rx=c \\Longrightarrow x``

这样针对对称正定矩阵的方程组求解，等于降低了一般的计算代价来实现。

# Example
```jldoctest
julia> A = Float64[4 -2 2;-2 2 -4; 2 -4 11]
3×3 Array{Float64,2}:
  4.0  -2.0   2.0
 -2.0   2.0  -4.0
  2.0  -4.0  11.0
julia> R = cholesky_decomposition(A, 3)
3×3 Array{Float64,2}:
 2.0  -1.0   1.0
 0.0   1.0  -3.0
 0.0   0.0   1.0
```
"""
function cholesky_decomposition(A, n)
    R = zeros(Float64, (n, n))
    for k = 1:n
        if A[k,k] < 0; return; end
        R[k,k] = sqrt(A[k,k])
        u = A[k+1:n, k]/R[k,k]
        u_t = transpose(u)
        R[k, k+1:n] = u_t
        A[k+1:n, k+1:n] = A[k+1:n, k+1:n] - u*u_t
    end
    return R
end
