"""
# 复合辛普森法则
    simpson_rule(f, a, b, m)

## 辛普森法则

辛普森法则使用通过三点的二次函数替换函数。

令``f(x)``是具有连续二阶导数的函数，定义在``[x_0, x_2]``上，使用拉格朗日公式，得到具有误差项的二阶插值多项式是:

```math
f(x) = y_0\\frac {(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)} + y_1 \\frac {(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)} + y_2 \\frac {(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)} \\\\ 
+ \\frac {(x-x_0)(x-x_1)(x-x_2)}{3!}f'''(c_x) = P(x) + E(x)
```

则在区间``[x_0, x_2]``上：
```math
\\int_{x_0}^{x_2}f(x)dx = \\int_{x_0}^{x_2}P(x)dx + \\int_{x_0}^{x_2} E(x)dx
```

令``h=x_2-x_1=x_1-x_0``

``\\qquad \\begin{aligned} \\int_{x_0}^{x_2}P(x)dx = y_0\\int_{x_0}^{x_2}\\frac {(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}dx + y_1 \\int_{x_0}^{x_2}\\frac {(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)}dx+ \\\\ 
y_2 \\int_{x_0}^{x_2} \\frac {(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}dx = \\frac h3 (y_0 + 4y_1 + y_2) \\end{aligned}``

如果``f^{(4)}(x)```存在并连续，那么误差为:
```math
\\int_{x_0}^{x_2} E(x)dx = -\\frac {h^5}{90}f^{(4)}(c)
```

**辛普森法则**：

```math
\\int_{x_0}^{x_2}f(x)dx = \\frac h3 (y_0 + 4y_1 + y_2)-\\frac {h^5}{90}f^{(4)}(c)
```

其中``h=x_2-x_0=x_1-x_0``,``c``在``x_0``和``x_2``之间

## 复合辛普森法则

同样辛普森法则局限在单一区间上操作，由于积分在区间的所有子区间上具有可加性，我们可以将区间分为多个小区间分别处理。但复合辛普森法则是按照每个长2h的子区间
```math
\\int_{x_2i}^{x_{2i+2}}f(x)dx = \\frac h3 (f(x_{2i}) + 4f(x_{2i+1}) + f(x_{2i+2})) - \\frac {h^5}{90}f^{(4)}(c_i)
```

**复合辛普森法则公式**:
    simpson_rule(f, a, b, m)
```math
\\int_{a}^{b}f(x)dx =  \\frac h3 (y_0 + y_{2m} + 4\\sum_{i=1}^m y_{2i+1} + 2\\sum_{i=1}^{m-1}y_{2i}) - \\frac {(b-a)h^4}{180}f^{(4)}(c)
```
其中 ``h=\\frac {(b-a)}{m}``,c在a和b之间

同样对于误差项根据中值定理，改写为(``m = \\frac {b-a}{2h}``)：
```math
\\sum_{i=0}^{m-1}\\frac {h^5}{90}f^{(4)}(c_i) = \\frac {h^5}{90}mf^{(4)}(c) = \\frac {(b-a)h^4}{180}f^{(4)}(c)
```

# Example
```jldoctest
julia> simpson_rule(x-> log(exp(1), x), 1, 2, 4)   # ln(x)  [1,2]   m=4    误差很小了
0.3862920434663129
```
"""
function simpson_rule(f, a, b, m)
    h = (b-a)/(2m)
    x_odd = a+h:2h:b-h
    x_egg = a+2h:2h:b-2h
    return (h/3)*(f(a)+f(b)+4*sum(f.(x_odd)) + 2*sum(f.(x_egg)))
end
