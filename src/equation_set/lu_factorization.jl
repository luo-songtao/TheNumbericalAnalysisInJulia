using LinearAlgebra

"""
# 基于LU分解表示朴素的高斯消元法
    lu_factorization(A, b, n)

将系数矩阵分解为下三角矩阵L和上三角矩阵U，再求解x，可降低运算次数。

LU分解使得在朴素的高斯消元法中的消去过程，向量b没有参与运算，直到回代过程才参与运算(两次)，相比而言，当n非常大时，计算次数会有明显的差异

注意：当前实现的算法如果遇到0主元会抛出异常并终止

```math
A = LU \\\\
Ax=b \\Longleftrightarrow LUx=b \\\\
Lc = b \\Rightarrow c \\\\
Ux=c \\Rightarrow x \\\\
```

# Arguments
- `A`: 表示系数矩阵A
- `b`: 表示常数项b
- `n`: 方程数

# Example
```jldoctest
julia> A = Float64[1 2 -1;2 1 -2;-3 1 1]
3×3 Array{Float64,2}:
  1.0  2.0  -1.0
  2.0  1.0  -2.0
 -3.0  1.0   1.0
julia> b = Float64[3; 3; -6]
3-element Array{Float64,1}:
  3.0
  3.0
 -6.0
julia> x = lu_factorization(A,b,3)
3×1 Array{Float64,2}:
 3.0
 1.0
 2.0
```
"""
function lu_factorization(A, b, n)
    # 1. LU分解过程
    for j = 1:n-1
        if abs(A[j,j]) < eps(1.0)    # 等价于2.0^-52
            error("Zero pivot encounted")
        end
        for i = j+1:n
            mult = A[i,j]/A[j,j]    # 行变换乘子
            for k = j+1:n
                A[i,k] = A[i,k] - mult * A[j,k]
            end
            A[i,j] = mult
        end
    end

    # 2. 回代步骤一 Lc = b 解出c
    for i = 1:n
        for j = 1:i-1
            b[i] = b[i] - A[i,j]*b[j]
        end
    end
    # b即是c

    x = zeros(n, 1)
    # 3. 回代步骤二 Ux=c，解出x
    for i = n:-1:1
        for j = i+1:n
            b[i] = b[i] - A[i,j]*x[j]
        end
        x[i] = b[i]/A[i,i]
    end
    return x
end

# 以下是未优化版本
# function lu_factorization(a, b, n)
#     L = Matrix{Float64}(I, n, n)
#     # 1. LU分解过程
#     for j = 1:n-1
#         if abs(a[j,j]) < eps(1.0)    # 等价于2.0^-52
#             error("Zero pivot encounted")
#         end
#         for i = j+1:n
#             mult = a[i,j]/a[j,j]    # 行变换乘子
#             for k = j+1:n
#                 a[i,k] = a[i,k] - mult * a[j,k]
#             end
#             a[i,j] = 0    # 将对应位置置为0，即将a矩阵变为U矩阵
#             L[i,j] = mult    # 将乘子放置如L矩阵对应位置
#         end
#     end

#     # 2. 回代步骤一 Lc = b 解出c
#     c = zeros(n, 1)
#     for i = 1:n
#         for j = 1:i-1
#             b[i] = b[i] - L[i,j]*c[j]
#         end
#         c[i] = b[i]/L[i,i]
#     end

#     U = a
#     x = zeros(n, 1)
#     # 3. 回代步骤二 Ux=c，解出x
#     for i = n:-1:1
#         for j = i+1:n
#             c[i] = c[i] - U[i,j]*x[j]
#         end
#         x[i] = c[i]/U[i,i]
#     end
#     return x
# end
