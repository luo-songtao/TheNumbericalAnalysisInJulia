"""
# 多变量牛顿方法
    mult_newton_method(f, df, x_0, k)

多变量牛顿方法是由单变量牛顿方法推演而来。将单变量情况下的导函数``f'``，定义为向量值函数``F``的**雅可比矩阵**。

``\\vec x=(u,v,w,···)``

``F(x) = 0``

即
``F(x) = \\left [ \\begin{matrix} f_1(x) \\\\ f_2(x) \\\\ f_3(x) \\\\ ··· \\end{matrix} \\right ] = \\left [ \\begin{matrix} f_1(u, v, w,···) \\\\ f_2(u, v, w,···) \\\\ f_3(u, v, w,···) \\\\ ··· \\end{matrix} \\right ] = 0``

## 雅可比矩阵DF(x)

``DF(x) = \\left [ \\begin{matrix} \\frac{\\partial f_1}{\\partial u} \\space \\frac{\\partial f_1}{\\partial v} \\space \\frac{\\partial f_1}{\\partial w} \\space · \\\\ \\frac{\\partial f_2}{\\partial u} \\space \\frac{\\partial f_2}{\\partial v} \\space \\frac{\\partial f_2}{\\partial w} \\space · \\\\ \\frac{\\partial f_3}{\\partial u} \\space \\frac{\\partial f_3}{\\partial v} \\space \\frac{\\partial f_3}{\\partial w} \\space · \\\\ ··· \\end{matrix} \\right ]``

令``x=r``是根，``x_0``是当前近似估计，根据泰勒一级展开，忽略二次余项

``\\begin{aligned} 0 = F(r) &\\approx F(x_0) + DF(x_0)\\cdot (r-x_0) \\\\ -DF(x_0)^{-1}F(x_0) &\\approx r-x_0 \\end{aligned}``

## 多变量牛顿方法迭代公式：
``\\begin{aligned} x_0 &= 初始估计 \\\\ x_{k+1} &= x_k- (DF(x_k))^{-1}F(x_k) \\space , \\space k=0,1,2,··· \\end{aligned}``

由于计算矩阵逆的代价通常较大，令``s``是``DF(x_k)s = F(x_k)``的解，然后在每步中使用高斯消去代替计算矩阵的逆，减少运算次数

``\\begin{aligned} &DF(x_k)s = F(x_k) \\\\ &x_{k+1} = x_k - s \\end{aligned}``

# Example
```jldoctest
# F = [(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)]
# DF = [(u,v)->(-3*u^2) (u,v)->(1); (u,v)->(2u) (u,v)->(2v)]
# x_0 = [1, 2]
julia> mult_newton_method([(u,v)->(v-u^3); (u,v)->(u^2 + v^2 -1)],  [(u,v)->(-3*u^2) (u,v)->(1); (u,v)->(2u) (u,v)->(2v)], [1, 2], 10)
2×1 Array{Float64,2}:
 0.826031357654187
 0.5636241621612585
```
"""
function mult_newton_method(F, DF, x_0, k)
    n = length(x_0)
    x_k = x_0
    for i = 1:k
        y_k = Float64[f(x_k...) for f in F]
        df_k = Float64[df(x_k...) for df in DF]
        reshape(df_k, (n, n))
        s = naive_gauss_elimination(df_k, y_k, n)
        x_k = x_k - s
    end
    return x_k
end

push!(LOAD_PATH,"../")
using EquationSet