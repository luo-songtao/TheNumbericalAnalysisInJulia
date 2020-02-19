"""
# 高斯-牛顿方法
    guass_newton(r, Dr, x_0, k)

对于非线性最小二乘问题，可以使用高斯牛顿方法求解。

## 高斯牛顿方法推导

目的是最小化误差余项的平方和，对于m个方程和n个未知数的方程组，误差为:

``\\qquad \\begin{aligned} r_1(x_1, x_2, ..., x_n) &= 0 \\\\ & ... \\\\ r_m(x_1, x_2, ..., x_n) &= 0 \\end{aligned}``

误差平方和为：

```math
E(x_1, ..., x_n) = \\frac 12 (r^2_1+...+r^2_m)
```

## 高斯牛顿方法推导源于牛顿法：

``E(x)``的二阶Taylor展开:

``\\qquad \\begin{aligned} E(x) \\approx E(x_k) + E'(x_k)(x-x_k) + \\frac 12 E''(x_k)(x-x_k)^2  \\end{aligned}``

令``\\Delta x = x - x_k``,即``x = x_k + \\Delta x``:

``\\qquad \\begin{aligned} E(x_k+\\Delta x) \\approx E(x_k) + E'(x_k)\\Delta x + \\frac 12 E''(x_k){\\Delta x}^2 \\end{aligned}``

为了最小化E，也就是极值问题，因此令E的梯度为0:

``\\qquad g(x_k) = \\nabla E(x_k) = \\frac {d}{d\\Delta x}(E(x_k) + E'(x_k)\\Delta x + \\frac 12 E''(x_k){\\Delta x}^2) = E'(x_k) + E''(x_k)\\Delta x = 0``

则有：

``\\Delta x = - \\frac {E'(x_k)}{E''(x_k)}``

由于E是向量函数，这里将``x_{k+1} = x_k + \\Delta x``改写为
```math
x_{k+1} = x_k - H^{-1}g
```
其中g为梯度向量，H为海森矩阵

``\\qquad \\begin{aligned} g_j &= 2\\sum^m_{i=1}r_i \\frac {\\partial r_i}{\\partial x_j} = 2\\sum^m_{i=1}r_i J_{ij} \\qquad & j=1,2,...,n \\\\ h_{jk} &= 2\\sum^m_{i=1}(\\frac{\\partial r_i}{\\partial x_j}\\frac{\\partial r_i}{\\partial x_k} + \\frac{\\partial^2 r_i}{\\partial x_j \\partial x_k}) \\approx 2\\sum^m_{i=1} J_{ij} J_{ik} \\qquad & j,k=1,2,...,n \\end{aligned}``

简化为:

``\\qquad \\begin{aligned} g &= 2J^T_rr \\\\ H &\\approx 2J^T_rJ_r \\end{aligned}``

其中``J_r``表示雅可比矩阵

```math
x_{k+1} = x_k - (J^T_rJ_r)^{-1}J^T_rr(x_k)
```

同样为了避免求逆的过程，令``H\\Delta{x_k} = -g``,即``J^T_rJ_r \\Delta{x_k} = -J^T_rr(x_k)``，先求解出``\\Delta{x_k}``,
然后代入``x_{k+1} = x_k + \\Delta{x_k}``更新

# Example
```julia
function test_guass_newton()
    # 找到一点，要求离三个圆的距离的平方和最小
    # 圆心坐标
    x1, y1 = -1, 0
    x2, y2 = 1, 0.5
    x3, y3 = 1, -0.5
    # 圆半径
    r1, r2, r3 = 1, 0.5, 0.5
    # 误差函数
    r = [
        (x,y) -> sqrt((x-x1)^2+ (y-y1)^2) - r1
        (x,y) -> sqrt((x-x2)^2+ (y-y2)^2) - r2
        (x,y) -> sqrt((x-x3)^2+ (y-y3)^2) - r3
    ]
    # r关于(x,y)的一阶偏导数矩阵
    Dr = [
        (x,y) -> (x-x1)/sqrt((x-x1)^2+ (y-y1)^2) (x,y) -> (y-y1)/sqrt((x-x1)^2+ (y-y1)^2);
        (x,y) -> (x-x2)/sqrt((x-x2)^2+ (y-y2)^2) (x,y) -> (y-y2)/sqrt((x-x2)^2+ (y-y2)^2);
        (x,y) -> (x-x3)/sqrt((x-x3)^2+ (y-y3)^2) (x,y) -> (y-y3)/sqrt((x-x3)^2+ (y-y3)^2);
    ]
    x = round.(guass_newton(r, Dr, [0;0], 10), digits=6)
    @assert reshape(x, 1,2) == [0.412891 0.0]
    return 
end
```
```jldoctest
julia> test_guass_newton()
2×1 Array{Float64,2}:
 0.412891
 0.0 
```
"""
function guass_newton(r, Dr, x_0, k)
    m,n = size(Dr)
    x_k = x_0
    for k = 1:k
        r_k = Float64[r_f(x_k...) for r_f in r]
        J_k = Float64[dr_f(x_k...) for dr_f in Dr]
        Δx = naive_gauss_elimination(J_k'J_k, -J_k'r_k, n)
        x_k = x_k + Δx
    end
    return x_k
end

push!(LOAD_PATH,"../")
using EquationSet

function test_guass_newton()
    # 找到一点，要求离三个圆的距离的平方和最小
    # 圆心坐标
    x1, y1 = -1, 0
    x2, y2 = 1, 0.5
    x3, y3 = 1, -0.5
    # 圆半径
    r1, r2, r3 = 1, 0.5, 0.5
    # 误差函数
    r = [
        (x,y) -> sqrt((x-x1)^2+ (y-y1)^2) - r1
        (x,y) -> sqrt((x-x2)^2+ (y-y2)^2) - r2
        (x,y) -> sqrt((x-x3)^2+ (y-y3)^2) - r3
    ]
    # r关于(x,y)的一阶偏导数矩阵
    Dr = [
        (x,y) -> (x-x1)/sqrt((x-x1)^2+ (y-y1)^2) (x,y) -> (y-y1)/sqrt((x-x1)^2+ (y-y1)^2);
        (x,y) -> (x-x2)/sqrt((x-x2)^2+ (y-y2)^2) (x,y) -> (y-y2)/sqrt((x-x2)^2+ (y-y2)^2);
        (x,y) -> (x-x3)/sqrt((x-x3)^2+ (y-y3)^2) (x,y) -> (y-y3)/sqrt((x-x3)^2+ (y-y3)^2);
    ]
    x = round.(guass_newton(r, Dr, [0;0], 10), digits=6)
    @assert reshape(x, 1,2) == [0.412891 0.0]
    return x
end