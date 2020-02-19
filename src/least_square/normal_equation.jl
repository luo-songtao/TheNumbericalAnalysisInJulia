"""
# 最小二乘与法线方程

## 最小二乘

- 求解方程时，有可能方程的个数超过未知变量个数或方程本身不存在解时，可以通过**最小二乘近似**找到第二可能好的解
- 在寻找多项式，并精确拟合数据点时，如果有大量的数据点，或者采集的数据点具有一定误差，使用高阶多项式精确拟合一般不是一个好方法。此时使用简单模型**近似拟合数**据是一种更合理的方式

在方程无解的情况下，可以通过一种直接方法找到最接近的``x``，这个特殊的``x``称为**最小二乘解**，我们这里把它记为: ``\\overline{x}``

也就是说``A\\overline{x} \\neq b``但最接近``b``，等价于：``b-A\\overline{x}``与平面``\\lbrace Ax|x\\in R^n \\rbrace``

## 法线方程

最小二乘基于**正交**。从平面外一点到一个平面的最短距离，有一个到平面的线段表示。**法线方程**可以确定该线段，它表示着最小二乘的误差。

已知：``(b-A\\overline{x})\\perp \\lbrace Ax|x\\in R^n \\rbrace``

将垂直性表示为矩阵的乘法，则可以发现对于``R^n``上的所有``x``:

``\\qquad \\begin{aligned} (Ax)^T(b-A\\overline{x}) &= 0 \\\\ x^TA^T(b-A\\overline{x}) &= 0 \\end{aligned}``

这也就意味着``n``维向量``A^T(b-A\\overline{x})``与中的其他维向量垂直，并且好包括自身，当且仅当:
```math
A^T(b-A\\overline{x}) = 0
```
即：
```math
A^TA\\overline{x} = A^Tb
```

这个方程组被称为**法线方程**，它的解``\\overline{x}``是方程组``Ax=b``的最小二乘解。它可以最小化余项``r=b-Ax``的欧式长度。

如果余项是0向量，那么意味着精确求解了``Ax=b``。否则余项向量的欧式长度是后向误差，它度量了``\\overline{x}``到解的距离

``\\qquad \\begin{aligned} \\Vert r \\Vert_2 &= \\sqrt{r^2_1 + ··· + r^2_m}，2范数 \\\\ SE &= r^2_1 + ··· + r^2_m，平方误差 \\\\ RMSE &= \\sqrt{\\frac{SE}{m}} = \\sqrt{r^2_1 + ··· + r^2_m}，均方根误差 \\end{aligned}``

## 使用法线方程求解最小二乘问题
```julia
A = Float64[1 -4; 2 3; 2 2]
b = Float64[-3;15;9]
x = round.(normal_equation(A, b), digits=2)
# [3.8; 1.8]
r = round.(b-A*x, digits=2)
# [0.4; 2.0; -2.2]
using LinearAlgebra
norm(r)
# 3.0
```
"""
function normal_equation(A, b)
    A_T = transpose(A)
    A = A_T*A
    b = A_T*b
    return lu_factorization_with_permutation(A, b, length(b))    # 使用PA=LU分解求解方程组
end

push!(LOAD_PATH, "../")
using EquationSet


"""
# 最小二乘线性建模

最小二乘建模：使用最小二乘对数据进行模型拟合

最小二乘核心思想：在数据点上通过平方误差度量拟合的余项，并找出模型的参数使得该误差最小

## 法线方程拟合数据

#### 拟合数据点：(-1, 1)、(0, 0)、(1, 0)、(2, -2)
```jldoctest
julia> x = [-1;0;1;2]
4-element Array{Int64,1}:
 -1
  0
  1
  2
julia> y = [1;0;0;-2]
4-element Array{Int64,1}:
  1
  0
  0
 -2
julia> fit_modeling_by_normal_equation([x.^0 x], x, y)    # 最优直线y = 0.2 - 0.9x
RMSE: 0.4183
2-element Array{Float64,1}:
  0.2
 -0.9
julia> fit_modeling_by_normal_equation([x.^0 x x.^2], x, y)    # 最优抛物线y=0.45 - 0.65x - 0.25x^2
RMSE: 0.3354
3-element Array{Float64,1}:
  0.45
 -0.65
 -0.25
```

## 周期数据拟合
```jldoctest
julia> t = 0:1/8:7/8
0.0:0.125:0.875

julia> y = [-2.2; -2.8; -6.1; -3.9; 0.0;1.1;-0.6;-1.1]
8-element Array{Float64,1}:
 -2.2
 -2.8
 -6.1
 -3.9
  0.0
  1.1
 -0.6
 -1.1

julia> fit_modeling_by_normal_equation([t.^0 cos.(t.*2pi) sin.(t.*2pi)], t, y)    # 模型： y = c1 + c2*cos2πt + c3*sin2πt
RMSE: 1.0631
3-element Array{Float64,1}:
 -1.95  
 -0.7445
 -2.5594

julia> fit_modeling_by_normal_equation([t.^0 cos.(t.*2pi) sin.(t.*2pi) cos.(t.*4pi) sin.(t.*4pi)], t, y)    # 模型： y = c1 + c2*cos2πt + c3*sin2πt + c4*cos4πt + c5*sin4πt
RMSE: 0.3962
5-element Array{Float64,1}:
 -1.95  
 -0.7445
 -2.5594
  1.125 
  0.825
```
"""
function fit_modeling_by_normal_equation(A, x, y)
    c = (A'*A)\(A'*y)
    r = y-A*c
    rmse = round(norm(r)/sqrt(length(r)), digits=4)
    println("RMSE: ", rmse)
    return round.(c, digits=4)
end
using LinearAlgebra


