"""
# 割线方法
    secant(f, x_0, x-1, k)
割线方法法计算方程解。
- 当难以计算函数导数时，割线方法是牛顿迭代法的一个很好的替代。
- 割线方法需要提供两个初始估计

## 割线方法迭代公式
``\\begin{aligned} x_0, x_1 &= 初始估计值\\\\ x_{i+1} &= x_i - \\frac {f(x_i)(x_i-x_{i-1})}{f(x_i)-f(x_{i-1})},(i=0,1,2,3,...) \\end{aligned}``

# Arguments
- `f`: `f(x)`函数
- `x_0`: 初始估计值1
- `x_1`: 初始估计值2
- `k`: 迭代次数

# Example
```jldoctest
julia> secant(x->x^3 +x -1, 0, 1, 20)
0.6823278038280193
```
"""
function secant(f, x_0, x_1, k)
    x_i_minus_1 = x_0
    x_i = x_1
    for i = 1:k
        # make sure the denominator not equals zero
        if f(x_i) == f(x_i_minus_1)
            break
        end
        x_i_plus_1 = x_i - f(x_i)*(x_i-x_i_minus_1)/(f(x_i)-f(x_i_minus_1))
        x_i_minus_1 = x_i
        x_i = x_i_plus_1
    end
    return x_i
end


