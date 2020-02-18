"""
# 格拉姆-施密特正交方法
    classic_gram_schmidt_orthogon(A, n)

格拉姆-施密特方法是对一组向量正交化。给定一组输入的m维向量，目的是找出正交坐标系统，获取由这些向量张成的空间。

更精确的说，给定n个线性无关的输入向量，该方法将计算出n个彼此正交的单位向量，构成一组单位正交向量。


## 格拉姆-施密特过程

令``a_1, a_2,···, a_n``是`R^m`中的线性无关向量``(n\\le m)``

首先定义：
```math
y_1=a_1 \\space 与 \\space q_1 = \\frac {y_1}{\\Vert y_1 \\Vert_2}
```
这样q_1是找出的第一个单位向量。为了找到第2个单位向量，对a_2在q_1方向上的投影减去a_2，并对结果规范化:

```math
y_2 = a_2 - q_1(q^T_1a_2), \\space q_2 = \\frac {y_2}{\\Vert y_2 \\Vert_2}
```
- tips：``(u,v) = v^Tu, \\space (u,u) = \\Vert u \\Vert^2``,``u、v``维数相同

则可以看到:

``\\qquad q^T_1y_2 = q^T_1a_2 - q^T_1a_2 = 0``

说明``q_1、q_2``是正交关系

因此定义第j步：
```math
y_j = a_j - q_1(q^T_1a_j) - q_2(q^T_2a_j)- ··· - q_{j-1}(q^T_{j-1}a_j), q_j = \\frac{y_j}{\\Vert y_j \\Vert_2}
```

这样的到的``q_1、q_2、···、q_n``两两正交，而且``q_1、q_2、···、q_n``张成的空间和``a_1, a_2,···, a_n``的等价     

# Example
```jldoctest
julia> A = [1 -4; 2 3; 2 2]
3×2 Array{Int64,2}:
 1  -4
 2   3
 2   2
julia> classic_gram_schmidt_orthogon(A, 2)
3×2 Array{Float64,2}:
 0.3333  -0.9333
 0.6667   0.3333
 0.6667   0.1333
```
"""
function classic_gram_schmidt_orthogon(A, n)
    q = zeros(length(A[:, 1]), n)
    for j = 1:n
        y = A[:, j]
        for i = 1:j-1
            r = q[:, i]'*A[:,j]
            y = y - r*q[:, i]
        end
        q[:, j] = y / sqrt(sum(y.^2))
    end
    return round.(q, digits=4)
end

"""
# 改进的格拉姆-施密特正交 
    gram_schmidt_orthogon(A, n)

对于格拉姆-施密特的微笑改进可以在及其计算中改进精度，改进的格拉姆-施密特正交与前面的原始方法在数学上等价

唯一不同的是内层循环中，``a_j``被``y``所替换。本质上它是等价的。

# Example
```jldoctest
julia> A = [1 1 1; 10^(-10) 0 0 ; 0 10^(-10) 0; 0 0 10^(-10)]
4×3 Array{Float64,2}:
 1.0      1.0      1.0    
 1.0e-10  0.0      0.0    
 0.0      1.0e-10  0.0    
 0.0      0.0      1.0e-10

julia> classic_gram_schmidt_orthogon(A, 3)    # 10^-20是一个可接受的双精度值，但在1 + 10^-20 = 1 ，此时精度丢失了
4×3 Array{Float64,2}:
 1.0   0.0      0.0   
 0.0  -0.7071  -0.7071
 0.0   0.7071   0.0   
 0.0   0.0      0.7071

julia> gram_schmidt_orthogon(A, 3)    # 这里一定程度避免了1 + 10^-20 这样的情况出现，从而改进了精度，得到准确的值
4×3 Array{Float64,2}:
 1.0   0.0      0.0   
 0.0  -0.7071  -0.4082
 0.0   0.7071  -0.4082
 0.0   0.0      0.8165
```
"""
function gram_schmidt_orthogon(A, n)
    q = zeros(Float64, length(A[:, 1]), n)
    for j = 1:n
        y = A[:, j]
        for i = 1:j-1
            r = q[:, i]'*y    # 修改的地方
            y = y - r*q[:, i]
        end
        q[:, j] = y / sqrt(sum(y.^2))
    end
    return round.(q, digits=4)
end
