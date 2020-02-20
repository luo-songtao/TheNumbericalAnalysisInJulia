"""
# 复合梯形法则
    trapezoidal_rule(f, a, b, m)

## 梯形法则

梯形法则使用通过``(x_0,f(x_0))和(x_1,f(x_1))``的直线替换函数。

令``f(x)``是具有连续二阶导数的函数，定义在``[x_0, x_1]``上，使用拉格朗日公式，得到具有误差项的一阶插值多项式是:

```math
f(x) = y_0\\frac {x-x_1}{x_0-x_1} + y_1 \\frac {x-x_0}{x_1-x_0} + \\frac {(x-x_0)(x-x_1)}{2!}f''(c_x) = P(x) + E(x)
```

则在区间``[x_0, x_1]``上：
```math
\\int_{x_0}^{x_1}f(x)dx = \\int_{x_0}^{x_1}P(x)dx + \\int_{x_0}^{x_1} E(x)dx
```

令``h=x_1-x_0``

``\\qquad \\begin{aligned} \\int_{x_0}^{x_1}P(x)dx &= y_0\\int_{x_0}^{x_1}\\frac {x-x_1}{x_0-x_1} + y_1 \\int_{x_0}^{x_1}\\frac {x-x_0}{x_1-x_0} = \\frac h2 (y_0+y_1) \\\\ \\int_{x_0}^{x_1}E(x)dx &= \\frac 1{2!} \\int_{x_0}^{x_1}(x-x_0)(x-x_1)f''(c)dx = -\\frac {h^3}{12}f''(c) \\end{aligned}``

**梯形法则公式**：

```math
\\int_{x_0}^{x_1}f(x)dx = \\frac h2 (y_0+y_1) - \\frac {h^3}{12}f''(c)
```

其中``h=x_1-x_0``,``c``在``x_0``和``x_1``之间

## 复合梯形法则

梯形法则局限在单一区间上操作，由于积分在区间的所有子区间上具有可加性，我们可以将区间分为多个小区间分别处理，这样的策略被称为复合数值积分

```math
\\int_{x_i}^{x_{i+1}}f(x)dx = \\frac h2 (f(x_i)+f(x_{i+1})) - \\frac {h^3}{12}f''(c_i)
```

**复合梯形法则公式**:

```math
\\int_{a}^{b}f(x)dx = \\frac h2 (y_0+y_m + 2\\sum_{i=1}^{m-1}y_i) - \\frac {(b-a)h^2}{12}f''(c)
```
其中 ``h=\\frac {(b-a)}{m}``,c在a和b之间

误差项根据中值定理，改写为：
```math
\\sum_{i=0}^{m-1}\\frac {h^3}{12}f''(c_i) = \\frac {h^3}{12}mf''(c) = \\frac {(b-a)h^2}{12}f''(c)
```

# Example
```jldoctest
julia> trapezoidal_rule(x-> log(exp(1), x), 1, 2, 1)   # ln(x)  [1,2]   m=1
0.34657359027997264
julia> trapezoidal_rule(x-> log(exp(1), x), 1, 2, 10)   # ln(x)  [1,2]  m=10
0.3216925481285146
```
"""
function trapezoidal_rule(f, a, b, m)
    h = (b-a)/m
    x = a:h:b
    if length(x[2:m-1]) == 0
        return (f(a) + f(b))*(h/2)
    else
        return (sum(2*f.(x[2:m-1])) + f(a) + f(b))*(h/2)
    end
end
