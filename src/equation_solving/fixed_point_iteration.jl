"""
# 不动点迭代法
    fixed_point_iteration(g, x_0, k)
不动点迭代法求解函数值

```math
f(x)=g(x)-x\\\\
x_0 = 初始估计\\\\
x_{i+1} = g(x_i),(i=0,1,2,3,...)
```

# Arguments
- `g`: `g(x)`函数
- `x_0`: 初始估计值
- `k`: 迭代次数

# Usage
```jldoctest
julia> fixed_point_iteration(x-> cos(x), 0, 100)
0.7390851332151607
```
"""
function fixed_point_iteration(g, x_0, k)
    x_i = x_0
    for i = 1:k
        x_i_plus_1 = g(x_i)
        x_i = x_i_plus_1
    end
    return x_i
end