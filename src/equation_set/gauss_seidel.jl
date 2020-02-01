"""
# 高斯-赛德尔方法
    gauss_seidel(A, b, x_0, n, k)

高斯-赛德尔方法和雅可比方法非常相似，唯一的差异的在于，高斯-赛德尔方法在每一步中都会用到最近更新的未知变量的值。

同样如果矩阵`A`是严格对角占优矩阵(主对角线的值在其所在行是最大的)，那么`A`是非奇异矩阵，且对于所有的向量`b`和初始估计，对`Ax=b`应用高斯-赛德尔方法都会收敛到(唯一)解

下式中：
- `D`表示`A`的主对角线矩阵
- `L`表示`A`的下三角矩阵(主对角线以下的元素)
- `U`表示`A`的上三角矩阵(主对角线以上的元素)

```math
Ax = b \\Longleftrightarrow (L+D)x_{k+1} = b - Ux_k
```

## 高斯-赛德尔方法迭代公式：
``x_0 =`` 初始估计
```math
x_{k+1} = D^{-1}(b-Lx_{k+1}-Ux_k), k=0,1,2,......
```

# Arguments
- `A`: 表示系数矩阵A
- `b`: 表示常数项b
- `x_0`: 初始估计(向量)
- `n`: 方程数
- `k`: 迭代次数

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
julia> x = gauss_seidel(A, b, x_0, 3, 50)    # 该例中约50步收敛
3-element Array{Float64,1}:
  2.0
 -1.0
  1.0
```
"""
function gauss_seidel(A, b, x_0, n, k)
    x_k = x_0
    for i = 1:k
        for j = 1:n
            L_mult = reshape(A[j, 1:j-1], (1, j-1)) * reshape(x_k[1:j-1], (j-1, 1))
            U_mult = reshape(A[j, j+1:n], (1, n-j)) * reshape(x_k[j+1:n], (n-j, 1))
            x_k[j] = (1/A[j,j]) * (b[j] - L_mult[1] - U_mult[1])
        end
    end
    return x_k
end
