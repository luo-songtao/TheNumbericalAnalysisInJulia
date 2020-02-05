"""
# 牛顿差商
    newton_difference_quotient
牛顿差商给出插值多项式的一种简单形式。给定``n``个数据点，所得到的结果多项式之多``n-1``阶

## 牛顿差商公式定义：

用``f[x_1 ··· x_n]``表示(唯一)多项式的下``x^{n-1}``项的系数，该多项式的插值``(x_1, f(x_1)),···,(x_n, f(x_n))``

``\\begin{aligned} P(x) = & f[x_1] + f[x_1 x_2](x-x_1) + f[x_1 x_2 x_3](x-x_1)(x-x_2) \\\\ &+ f[x_1 x_2 x_3 x_4](x-x_1)(x-x_2)(x-x_3) \\\\ &+ ··· + f[x_1 ··· x_n](x-x_1)···(x-x_{n-1}) \\end{aligned}``

根据唯一性，``x_1 x_2 ··· x_n``的任意置换结果相同：

``f[x_1 x_2 ··· x_n] = f[x_2 x_3 ··· x_n x_1] = f[x_2 x_3 ···x_{n-1} x_1 x_n]``

由此可推导出:

``f[x_1 ··· x_k] = \\frac {f[x_2 ··· x_k] - f[x_1 ··· x_{k-1}]}{x_k - x_1}``

以及：

``\\begin{aligned} f[x_k] &= f(x_k) \\\\ f[x_k x_{k+1}] &= \\frac {f[x_{k+1}]-f[x_k]}{x_{k+1}- x_{k}} \\\\ f[x_k x_{k+1} ··· x_{k+j}] &= \\frac {f[x_{k+1} ··· x_{k+j}]-f[x_k ··· x_{k+j-1}]}{x_{k+j}- x_{k}} \\end{aligned}``

# Example
```jldoctest
julia> x = [0 2 3]
1×3 Array{Int64,2}:
 0  2  3
julia> y = [1 2 4]
1×3 Array{Int64,2}:
 1  2  4
julia> newton_difference_quotient(x,y)
3-element Array{Float64,1}:
 1.0
 0.5
 0.5
```

```julia
using Plots
x_0 = 0:0.01:4
y_0 = horner_rule(2, c, x_0, x)
pyplot()
p = scatter(x,y)
plot!(p, x_0,y_0)
```
![](../img/newton_difference_quotient.png)

"""
function newton_difference_quotient(x, y)
    n = length(x)
    # coefficients[i,j]表示f[x_i ··· x_j]
    coefficients = zeros(n,n)

    for i = 1:n
        coefficients[i,i] = y[i]
    end

    for l = 2:n
        for i = 1:n-l+1
            coefficients[i,i+l-1] = (coefficients[i+1,i+l-1] - coefficients[i,i+l-2])/(x[i+l-1]-x[i])
        end
    end
    # 最后第一行将分别是f[x_1]、f[x_1 x_2]、···、f[x_1 ··· x_n]的值
    return coefficients[1, :]
end

x = [0 2 3]
y = [1 2 4]

c = newton_difference_quotient(x,y)


