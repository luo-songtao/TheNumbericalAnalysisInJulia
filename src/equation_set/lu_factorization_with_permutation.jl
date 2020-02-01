using LinearAlgebra

"""
# 使用置换的LU分解
    lu_factorization_with_permutation(a, b, n)

为避免淹没问题，在对矩阵A进行LU分解前，先对主元列进行判断，将主元列最大的一行与当前的首行进行置换

```math
P: 置换矩阵 \\\\
A = LU \\Longleftrightarrow PA=LU \\\\
```

# Arguments
- `a`: 表示系数矩阵A
- `b`: 表示常数项b
- `n`: 方程数

# Usage
```jldoctest
julia> m = Float64[2 1 5; 4 4 -4; 1 3 1]
3×3 Array{Float64,2}:
 2.0  1.0   5.0
 4.0  4.0  -4.0
 1.0  3.0   1.0
julia> b = Float64[5; 0; 6]
3-element Array{Float64,1}:
 5.0
 0.0
 6.0
julia> x = lu_factorization_with_permutation(m,b,3)
3×1 Array{Float64,2}:
 -1.0
  2.0
  1.0
```
"""
function lu_factorization_with_permutation(a, b, n)
    # 1. LU分解过程
    for j = 1:n-1
        if abs(a[j,j]) < eps(1.0)    # 等价于2.0^-52
            error("Zero pivot encounted")
        end
        # 向后找出当前列中的最大行
        max_row = j
        for temp = j+1:n
            if a[temp,j] > a[max_row, j]
                max_row = temp
            end
        end
        # 进行行置换
        if max_row != j
            a[j, :], a[max_row, :] = a[max_row, :], a[j, :]
            b[j], b[max_row] = b[max_row], b[j]
        end

        for i = j+1:n
            mult = a[i,j]/a[j,j]    # 行变换乘子
            for k = j+1:n
                a[i,k] = a[i,k] - mult * a[j,k]
            end
            a[i,j] = mult
        end
    end

    # 2. 回代步骤一 Lc = b 解出c
    for i = 1:n
        for j = 1:i-1
            b[i] = b[i] - a[i,j]*b[j]
        end
    end
    # b即是c

    x = zeros(n, 1)
    # 3. 回代步骤二 Ux=c，解出x
    for i = n:-1:1
        for j = i+1:n
            b[i] = b[i] - a[i,j]*x[j]
        end
        x[i] = b[i]/a[i,i]
    end
    return x
end

