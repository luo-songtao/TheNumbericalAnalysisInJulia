"""
# 雅可比方法
    jacobi(a, b, x_0, n, k)

雅可比方法是方程组系统中的一种形式的不动点迭代。在FPI中第一步是重写方程，进而求解未知量。

雅可比方法中，是如下步骤：
- 求解第i个方程得到第i个未知变量
- 然后使用不动点迭代，从初始估计开始，进行迭代

如果矩阵`A`是严格对角占优矩阵(主对角线的值在其所在行是最大的)，那么`A`是非奇异矩阵，且对于所有的向量`b`和初始估计，对`Ax=b`应用雅可比方法都会收敛到(唯一)解

下式中：
- `D`表示`A`的主对角线矩阵
- `L`表示`A`的下三角矩阵(主对角线以下的元素)
- `U`表示`A`的上三角矩阵(主对角线以上的元素)

```math
Ax = b \\Longleftrightarrow (L+D+U)x = b \\\\
Dx = b-(L+U)x \\\\
x = D^{-1}(b-(L+U)x)
```
## 雅可比方法迭代公式：
``x_0 =`` 初始估计
```math
x_{k+1} = D^{-1}(b-(L+U)x_k), k=0,1,2,......
```

# Arguments
- `a`: 表示系数矩阵A
- `b`: 表示常数项b
- `x_0`: 初始估计(向量)
- `n`: 方程数
- `k`: 迭代次数

# Usage
```jldoctest
julia> m = Float64[3 1; 1 2]
2×2 Array{Float64,2}:
 3.0  1.0
 1.0  2.0
julia> b = Float64[5; 5]
2-element Array{Float64,1}:
 5.0
 5.0
julia> x = jacobi(m,b,[0;0], 2, 100)
2-element Array{Float64,1}:
 1.0
 2.0
```
"""
function jacobi(a, b, x_0, n, k)
    D_inverse = zeros(Float64, (n, n))
    for i = 1:n
        D_inverse[i,i] = 1/a[i,i]    # D是对角线矩阵，其逆矩阵直接取倒
    end
    x_k = x_0
    for i = 1:k
        x_k = D_inverse * (b - (a-inv(D_inverse))*x_k)
    end
    return x_k
end
