"""
# 自适应积分
虽然对于复合的梯形法则和复合的辛普森法则等积分近似方法，可以得到对应的误差公式，但是这些误差公式都是依赖于更高阶的导数进行求取。
因此如果要想直接使用这些公式计算，用于满足给定容差对应的h值通常是比较困难的。

另外这些方法都是使用的相同的步长，，而函数通常啃根在定义域的某些部分的变化非常剧烈，而在其他部分变化缓慢，因此相同的步长明显是太适合整个区间。

因此这里使用积分误差公式，可以在运算中推出一个标准，其步长对于特定的子空间适合，这种方法称为自适应积分。

## 自适应梯形法则公式
    trap_adaptive_integral(f, a, b, tol)

在区间[a,b]上的梯形法则``S_{[a,b]}``满足：

```math
\\int_a^bf(x)dx = S_{[a,b]} - h^3\\frac {f''(c_0)}{12}
```

我们对区间[a,b]进行二分为[a,c],[c,b]后(c为中点)可得：

```math
\\int_a^bf(x)dx = S_{[a,c]} - \\frac {h^3}{8} \\frac {f''(c_1)}{12} + S_{[c,b]} - \\frac {h^3}{8} \\frac {f''(c_2)}{12} = S_{[a,c]} + S_{[c,b]} - \\frac {h^3}{4} \\frac {f''(c_3)}{12}
```

将上面的两式相减得到:

```math
S_{[a,b]} - (S_{[a,c]} + S_{[c,b]}) = - \\frac {h^3}{4} \\frac {f''(c_3)}{12} + h^3\\frac {f''(c_0)}{12} \\approx \\frac {3}{4}h^3 \\frac {f''(c_0)}{12}
```
这里我们近似``f''(c_0) \\approx f''(c_3)``

因此可以看到``S_{[a,b]} - (S_{[a,c]} + S_{[c,b]})``的误差，近似是``S_{[a,c]} + S_{[c,b]}``的误差的三倍

因此那么可以通过检查``S_{[a,b]} - (S_{[a,c]} + S_{[c,b]})``的值是小于给定容差的三倍，来作为一种近似方式，判断误差是否近似低于给定容差

如果没有满足，那么就对剩下的区间再进行递归形式划分，然后求解
# Example
使用自适应梯形积分近似计算``\\int_{-1}^1 (1+\\sin e^{3x})dx``
```jldoctest
julia> trap_adaptive_integral(x -> 1 + sin(exp(1)^(3x)), -1, 1, 0.005)
2.501859632312862
julia> trap_adaptive_integral(x -> 1 + sin(exp(1)^(3x)), -1, 1, 0.0005)
2.5009058956107495
```
"""
function trap_adaptive_integral(f, a, b, tol)
    c = (a+b)/2
    s_ab = trapezoidal_rule(f, a, b, 1)
    s_ac = trapezoidal_rule(f, a, c, 1)
    s_cb = trapezoidal_rule(f, c, b, 1)
    if abs(s_ab - (s_ac+s_cb)) < 3*tol
        return s_ac + s_cb
    else
        return trap_adaptive_integral(f, a, c, tol/2) + trap_adaptive_integral(f, c, b, tol/2)
    end
end

push!(LOAD_PATH,"../")
using NumericalCalculus


"""
## 自适应辛普森法则公式
    simpson_adaptive_integral(f, a, b, tol)
在区间[a,b]上的梯形法则``S_{[a,b]}``满足：

```math
\\int_a^bf(x)dx = S_{[a,b]} - \\frac {h^5}{90} f^{(4)}(c_0)
```

我们对区间[a,b]进行二分为[a,c],[c,b]后(c为中点)可得：

```math
\\int_a^bf(x)dx = S_{[a,c]} - \\frac {h^5}{32} \\frac {f^{(4)}(c_1)}{90} + S_{[c,b]} - \\frac {h^5}{32} \\frac {f^{(4)}(c_2)}{90} = S_{[a,c]} + S_{[c,b]} - \\frac {h^5}{16} \\frac {f^{(4)}(c_3)}{90}
```

将上面的两式相减得到:

```math
S_{[a,b]} - (S_{[a,c]} + S_{[c,b]}) = \\frac {h^5}{90} f^{(4)}(c_0) - \\frac {h^5}{16} \\frac {f^{(4)}(c_3)}{90} \\approx \\frac {15}{16} h^5 \\frac {f^{(4)}(c_3)}{90}
```
这里我们近似``f^{(4)}(c_0) \\approx f^{(4)}(c_3)``

因此可以看到``S_{[a,b]} - (S_{[a,c]} + S_{[c,b]})``的误差，近似是``S_{[a,c]} + S_{[c,b]}``的误差的15倍，但通常，将15设置为10，显得更保守点
# Example
使用自适应辛普森积分近似计算``\\int_{-1}^1 (1+\\sin e^{3x})dx``
```jldoctest
julia> simpson_adaptive_integral(x -> 1 + sin(exp(1)^(3x)), -1, 1, 0.005)
2.5002123674678263
julia> simpson_adaptive_integral(x -> 1 + sin(exp(1)^(3x)), -1, 1, 0.0005)
2.499918700559637
```
"""
function simpson_adaptive_integral(f, a, b, tol)
    c = (a+b)/2
    s_ab = simpson_rule(f, a, b, 1)
    s_ac = simpson_rule(f, a, c, 1)
    s_cb = simpson_rule(f, c, b, 1)
    if abs(s_ab - (s_ac+s_cb)) < 15*tol
        return s_ac + s_cb
    else
        return simpson_adaptive_integral(f, a, c, tol/2) + simpson_adaptive_integral(f, c, b, tol/2)
    end
end

push!(LOAD_PATH,"../")
using NumericalCalculus

