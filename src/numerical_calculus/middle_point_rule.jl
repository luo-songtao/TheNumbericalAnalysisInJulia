"""
# 复合中点法则

梯形法则和辛普森法则都被称为闭牛顿-科特斯方法，因为它们都是在闭区间上的，因此如果是开区间上，可以使用中点法则

中点法则其实只用了一个点来进行插值，且它相当于一条平行于x轴的直线，因此它的一阶导数为0

同样对于函数f，二阶导数在给定区间上连续。那么f在区间中点``\\omega = x_0+h/2，h=x_1-x_0``的一阶Taylor展开：
```math
f(x) = f(\\omega) + f(x-\\omega)f'(\\omega) +\\frac 12 (x-\\omega)^2 f''(c_x)
```

队上式两则进行积分可得**中点法则公式**：
```math
\\int_{x_0}^{x_1} f(x)dx = (x_1-x_0)f(\\omega) + 0 + \\frac 12 \\int_{x_0}^{x_1}f''(c_x)(x-\\omega)^2dx = hf(\\omega) + \\frac {h^3}{24}f''(c)
```

**复合中点法则公式**
```math
\\int_{a}^{b}f(x)dx = h\\sum_{i=1}^{m}f(\\omega_i) + \\frac {b-a}{24}h^2f''(c)
```
其中``h=\\frac {b-a}{m}``

# Example
```jldoctest
julia> middle_point_rule(x->sin(x)/x, 0, 1, 10)
0.9462085788431455
```
"""
function middle_point_rule(f, a, b, m)
    h = (b-a)/m
    ω = a+h/2:h:b
    return sum(h*f.(ω))
end
