"""
# 豪斯霍尔德反射子
    householder_reflector(x, ω)
尽管改进的格拉姆-施密特正交是计算矩阵的QR分解的有效方式，但它还是不是最好的方式。

豪斯霍尔德反射方法是更好的一种方法，它需要更少的计算，同时在舍入误差放大的意义上将也更稳定。

豪斯霍尔德反射子是一个正交矩阵，它通过``m-1``维平面反射``m``维向量

假设``x``和``\\omega``是具有相同欧几里得长度的向量，``\\Vert x \\Vert_2 = \\Vert \\omega \\Vert_2``,则``\\omega - x ``和``\\omega + x ``正交

定义向量``v = \\omega - x ``，考虑投影矩阵：
```math
P = \\frac {vv^T}{v^Tv}
```

投影矩阵性质: 

``\\qquad P^ - 2 = P``

``\\qquad Pv = v``

令``H = I - 2P``，则根据``v = \\omega - x ``和``\\omega + x ``正交，可证

```math
Hx = Ix - 2Px = \\omega - v -2\\frac {vv^T}{v^Tv}x = \\omega - \\frac {vv^T}{v^Tv}(\\omega + x) = \\omega
```

矩阵``H``被称为**豪斯霍尔德反射子**，且``H``是对称并正交的矩阵:

``\\qquad H^TH = HH = (I-2P)(I-2P)=I-4P+4P^2 = I``

# Example
```jldoctest
julia> x = [3;4]
2-element Array{Int64,1}:
 3
 4

julia> ω = [5;0]
2-element Array{Int64,1}:
 5
 0

julia> householder_reflector(x,ω)
2×2 Array{Float64,2}:
 0.6   0.8
 0.8  -0.6
```
"""
function householder_reflector(x, ω)
    n = length(x)
    v = ω - x
    P = v*v'/(v'v)
    return Matrix(I, n, n) - 2P
end
using LinearAlgebra


"""
## 豪斯霍尔德反射子与QR分解
    qr_by_householder_reflector(A)
豪斯霍尔德反射子可以用于推导QR分解的一个新的方法。而且该方法是典型矩阵进行QR分解的常用方法。

思想：将列向量x移动到坐标轴，并以此将0放在矩阵中。

一开始令``x_1``是``A``的第一列，令``\\omega = \\pm (\\Vert x_1 \\Vert_2, 0,0,···,0)``为第一个坐标轴上的向量，它们的欧几里得擦和高难度相同（理论上那种符号都可以，但为了数值稳定，一般将``x``的第一个元素选为正号，以避免两个近似相等的数字相减）

然后生成豪斯霍尔德反射子``H_1``并满足``H_1x= \\omega``。

这是对于``H_1A``矩阵的第一列除第一个元素外都是0，现在继续这种方式生成``H_2、H_3、···```,直到将``A``变为上三角，这样也就得到了``R``矩阵。同时有(对于4x3的矩阵A)：

```math
H_3H_2H_1A = R
```
由于豪斯霍尔德反射子是对称正交矩阵，则有：
```math
A = H_1H_2H_3R 
```
也就是说``Q = H_1H_2H_3``,则实现了``A=QR``分解

#Example
```jldoctest
julia> A = [1 -4; 2 3; 2 2] 
3×2 Array{Int64,2}:
 1  -4
 2   3
 2   2
julia> q,r = qr_by_householder_reflector(A)
([0.33333333333333337 -0.9333333333333336 -0.1333333333333333; 0.6666666666666666 0.3333333333333332 -0.6666666666666667; 0.6666666666666666 0.13333333333333347 0.7333333333333334], [3.0 1.9999999999999998; 0.0 4.999999999999999; 0.0 0.0])
julia> q
3×3 Array{Float64,2}:
 0.333333  -0.933333  -0.133333
 0.666667   0.333333  -0.666667
 0.666667   0.133333   0.733333
julia> r
3×2 Array{Float64,2}:
 3.0  2.0
 0.0  5.0
 0.0  0.0
```
"""
function qr_by_householder_reflector(A)
    m,n = size(A)
    R = zeros(m,n)
    Q = I
    reflector = I
    for j = 1:n
        x_j = reflector*A[:, j]
        ω_j = zeros(m-j+1)
        ω_j[1] = sqrt(sum(x_j[j:m].^2))
        reflector = householder_reflector(x_j[j:m], ω_j)
        R[1:j-1, j] = x_j[1:j-1]
        R[j:m, j] = ω_j
        H = hcat(vcat(I, zeros(m-j+1, j-1)), vcat(zeros(j-1,m-j+1), reflector))
        Q = Q*H
        if j == n
            return Q, R
        end
    end
    return Q, R
end

