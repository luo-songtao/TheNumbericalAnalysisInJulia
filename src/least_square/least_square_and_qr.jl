"""
# 不完全的QR分解
    incomplete_qr(A, m, n)

将格拉姆-施密特的正交过程写成矩阵形式，也就是将得到的所有单位向量作为矩阵Q的列向量，将每次正交的余项放入R矩阵右上角，具体如下：

令``r_{jj} = \\Vert y_j \\Vert_2``, ``r_{ij} = q^T_ia_j``

``\\qquad \\begin{aligned} a_1 &= r_{11}q_1 \\\\ a_2 &= r_{12}q_1 + r_{22}q_2 \\end{aligned}``

更一般地：

``\\qquad a_j = r_{1j}q_1 + ··· + r_{j-1,j}q_{j-1} + r_{jj}q_j``

```math
A = QR
```

``A``是``m\\times n``矩阵, ``Q``是``m\\times n``矩阵, ``R``是``n\\times n``上三角方阵

因为只有n个m维单位向量，无法张成出``R^m``空间，所以称之为不完全的QR分解

# Example

```julia
julia> A = [1 -4; 2 3; 2 2]
3×2 Array{Int64,2}:
 1  -4
 2   3
 2   2
julia> q,r = incomplete_qr(A, 3, 2)
3×2 Array{Float64,2}:
 0.3333  -0.9333
 0.6667   0.3333
 0.6667   0.1333
2×2 Array{Float64,2}:
 3.0  2.0
 0.0  5.0
```

"""
function incomplete_qr(A, m, n)
    q = zeros(Float64, m, n)
    r = zeros(n, n)
    for j = 1:n
        y = A[:, j]
        for i = 1:j-1
            r[i,j] = q[:, i]'*y
            y = y - r[i,j]*q[:, i]
        end
        r[j,j] = sqrt(sum(y.^2))
        q[:, j] = y / r[j,j]
    end
    return round.(q, digits=4), round.(r, digits=4)
end

"""
# 完全QR分解
    complete_qr(A, m, n)
```math
A = QR
```

``A``是``m\\times n``矩阵, ``Q``是``m\\times m``正交方阵, ``R``是``m\\times n``上三角矩阵

为了实现完全QR分解，当``n<m``时，需要在``A``中加上额外的``m-n``个向量(也是线性无关的)，从而使得``Q``可以张成``R^m``

**注意：使用格拉姆-施密特方法对``m\\times m``矩阵做QR分解，计算次数比LU分解的三倍还多，此外还有大约相同数量的加法**

# Example
```julia
julia> A = [1 -4; 2 3; 2 2]
3×2 Array{Int64,2}:
 1  -4
 2   3
 2   2
julia> A = hcat(A, [1;0;0])    # 增加额外的线性无关向量
3×3 Array{Int64,2}:
 1  -4  1
 2   3  0
 2   2  0
julia> q,r = complete_qr(A, 3, 2)
3×3 Array{Float64,2}:
 0.3333  -0.9333   0.1333
 0.6667   0.3333   0.6667
 0.6667   0.1333  -0.7333
3×2 Array{Float64,2}:
 3.0  2.0
 0.0  5.0
 0.0  0.0
```
"""
function complete_qr(A, m, n)
    q = zeros(Float64, m, m)
    r = zeros(m, n)
    for j = 1:n
        y = A[:, j]
        for i = 1:j-1
            r[i,j] = q[:, i]'*y
            y = y - r[i,j]*q[:, i]
        end
        r[j,j] = sqrt(sum(y.^2))
        q[:, j] = y / r[j,j]
    end

    for j = n+1:m
        y = A[:, j]
        for i = 1:j-1
            y = y - (q[:, i]'*y)*q[:, i]
        end
        q[:, j] = y / sqrt(sum(y.^2))
    end

    return round.(q, digits=4), round.(r, digits=4)
end

"""
# 通过QR分解实现最小二乘
    least_square_by_complete_qr(A, b, x_0)

A是``m\\times n``矩阵，求解``x``:
```math
Ax = b
```

找出完全QR分解：``A=QR``

令：
- ``\\widehat{R}= R``的上``n\\times n``子矩阵
- ``\\widehat{d} = d = Q^Tb``的上``n``个元素

求解``R\\overline{x} = Q^Tb``得到最小二乘解``\\overline{x}``

# Example
```julia
julia> A = [1 -4; 2 3; 2 2]
3×2 Array{Int64,2}:
 1  -4
 2   3
 2   2

julia> A = hcat(A, [1;0;0])    # 增加额外的向量
3×3 Array{Int64,2}:
 1  -4  1
 2   3  0
 2   2  0

julia> b = [-3;15;9]
3-element Array{Int64,1}:
 -3
 15
  9

julia> least_square_by_complete_qr(A, b, [1.0;1.0])
2-element Array{Float64,1}:
 3.8
 1.8
```
"""
function least_square_by_complete_qr(A, b, x_0)
    m = length(b)
    n = length(x_0)
    q, r = complete_qr(A, m, n)
    r_τ = r[1:n, :]
    d_τ = (q'*b)[1:n]
    x = x_0
    for i = n:-1:1
        for j = i+1:n
            d_τ[i] = d_τ[i] - r_τ[i,j]*x[j]
        end
        x[i] = d_τ[i]/r_τ[i,i]
    end
    return round.(x, digits=2)
end
