"""
# 自然三次样条
    natural_cubic_spline(x, y, k=10)

## 三次样条

给定``n``个点``(x_1, y_1), (x_2, y_2),···,(x_n, y_n)``,其中``x_i``不同，并且升序。通过点``(x_1, y_1), (x_2, y_2),···,(x_n, y_n)``的**三次样条**S(x)是一组三次多项式：

``\\qquad \\begin{aligned} S_1(x) &= y_1 + b_1(x-x_1) + c_1(x-x_1)^2 + d_1(x-x_1)^3 \\space\\space 在区间[x_1, x_2]上 \\\\ S_2(x) &= y_2 + b_2(x-x_2) + c_2(x-x_2)^2 + d_2(x-x_2)^3 \\space\\space 在区间[x_2, x_3]上 \\\\ &··· \\\\ S_{n-1}(x) &= y_{n-1} + b_{n-1}(x-x_{n-1}) + c_{n-1}(x-x_{n-1})^2 + d_{n-1}(x-x_{n-1})^3 \\space\\space 在区间[x_{n-1}, x_n]上 \\end{aligned}``

且具有以下性质：
- **性质1**：``S_i(x_i) = y_i,S_{i+1}(x_{i+1}) = y_{i+1},i=1,2,···,n-1```
- **性质2**: ``S'_{i-1}(x_i) = S'_{i}(x_i),i=2,···,n-1```
- **性质3**: ``S''_{i-1}(x_i) = S''_{i}(x_i),i=2,···,n-1```

意即相邻的曲线段在节点上具有相同的1阶、2阶导数

但注意在每一段中，满足这三个性质的三次多项式是无穷多个的。通常需要添加额外的约束条件来将问题域约束到唯一解。

根据性质1,可得出以下``n-1``个方程:

``\\qquad \\begin{aligned} y_2 &= S_1(x_2) = y_1 + b_1(x_2-x_1) + c_1(x_2-x_1)^2 + d_1(x_2-x_1)^3 \\\\ & ··· \\\\ y_n &= S_{n-1}(x_n) = y_{n-1} + b_{n-1}(x_n-x_{n-1}) + c_{n-1}(x_n-x_{n-1})^2 + d_{n-1}(x_n-x_{n-1})^3 \\end{aligned}``

根据性质2，可得出以下``n-2``个方程:

``\\qquad \\begin{aligned} 0 &= S'_1(x_2) - S'_2(x_2) = b_1 + 2c_1(x_2-x_1)+3d_1(x_2-x_1)^2 - b_2 \\\\ &··· \\\\ 0 &= S'_{n-2}(x_{n-1}) - S'_{n-1}(x_{n-1}) = b_{n-2} + 2c_{n-2}(x_{n-1}-x_{n-2})+3d_{n-2}(x_{n-1}-x_{n-2})^2 - b_{n-1} \\end{aligned}``

根据性质3，可得出以下``n-2``个方程:

``\\qquad \\begin{aligned} 0 &= S''_1(x_2) - S''_2(x_2) = 2c_1 + 6d_1(x_2-x_1) - 2c_2 \\\\ &··· \\\\ 0 &= S''_{n-2}(x_{n-1}) - S''_{n-1}(x_{n-1}) = 2c_{n-2} + 6d_{n-2}(x_{n-1}-x_{n-2}) - 2c_{n-1} \\end{aligned}``

共计``3n-5``个方程、``3n-3``个系数(b_i、c_i、d_i)

为了简化方程，引入额外的未知变量``c_n = S''_{n-1}(x_n)/2``,会使计算更简单。同时引入：``\\delta_i = x_{i+1}-x_i``、``\\Delta_i = y_{i+1}-y_i``

则根据性质3的公式：

```math
d_i = \\frac{c_{i+1}-c_i}{3\\delta_i} ,i=1,2,···,n-1
```

带入性质1的公式：

```math
b_i = \\frac{\\Delta_i}{\\delta_i} - c_i\\delta_i - d_i\\delta^2_i = \\frac{\\Delta_i}{\\delta_i} - \\frac{\\delta_i}{3}(2c_i+c_{i+1}) i=1,2,···,n-1
```

将以上两个公式带入性质2的公式中,则简化处以下``n-2``个方程：

``\\qquad \\begin{aligned} \\delta_1c_1 + 2(\\delta_1+\\delta_2)c_2 + \\delta_2c_3 &= 3(\\frac{\\Delta_2}{\\delta_2} - \\frac{\\Delta_1}{\\delta_1}) \\\\ &··· \\\\ \\delta_{n-2}c_{n-2} + 2(\\delta_{n-2}+\\delta_{n-1})c_{n-1} + \\delta_{n-1}c_n &= 3(\\frac{\\Delta_{n-1}}{\\delta_{n-1}} - \\frac{\\Delta_{n-2}}{\\delta_{n-2}}) \\end{aligned}``

## 自然三次样条

自然三次样条也就是在性质1~3的基础上添加以下性质：
- ``S''_1(x_1)=0``
- ``S''_{n-1}(x_n)=0``

将样条的开始和结束端点设定为拐点。类似这种在区间两个端点添加的附加条件，被称为**边界条件**

自然三次样条的条件可以得到另外两个方程：

```math
S''_1(x_1)) = 2c_1 = 0 \\\\
S''_{n-1}(x_n) = c_n = 0
```

结合三次样条关于``c_i``的``n-2``个方程组合为一个针对n个未知量``(c_i)``的n个方程，写成矩阵形式，可以看到它是一个**严格对角占优矩阵**，通过它解出``[c_1 ··· c_n]``，然后带入分别求解出``b_i、d_i``

# Example
```jldoctest
julia> natural_cubic_spline([0 1 2], [3 -2 1], 10)
3×3 Array{Float64,2}:
  2.0  0.0  -7.0
 -2.0  6.0  -1.0
  0.0  0.0   0.0
```
"""
function natural_cubic_spline(x, y, k=10)
    n = length(x)
    # Ac = r
    A = zeros(n,n)
    c_0 = zeros(n)
    r = zeros(n)

    δ = zeros(n-1)
    Δ = zeros(n-1)

    for i = 1:n-1
        δ[i] = x[i+1] - x[i]
        Δ[i] = y[i+1] - y[i]
    end

    for i = 2:n-1
        A[i, i] = 2(δ[i-1] + δ[i])
        A[i, i-1] = δ[i-1]
        A[i, i+1] = δ[i]
        r[i] = 3(Δ[i]/δ[i] - Δ[i-1]/δ[i-1])
    end

    A[1, 1] = 1
    A[n, n] = 1
    c = gauss_seidel(A, r, c_0, n, k)    # 使用高斯-塞德尔方法求解c

    b = zeros(n)
    d = zeros(n)
    for i = 1:n-1
        b[i] = (c[i+1] - c[i]) / (3δ[i])
        d[i] = Δ[i]/δ[i] - (δ[i]/3)*(2c[i] + c[i+1])
    end
    return [b c d]
end

push!(LOAD_PATH,"../")
using EquationSet
