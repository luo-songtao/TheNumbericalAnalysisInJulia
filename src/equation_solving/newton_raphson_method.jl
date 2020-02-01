"""
# 牛顿拉普森方法
    newton_raphson_method(f, fd, x_0, k, m=1)
牛顿拉普森迭代法计算方程解

```math
x_0 = 初始估计\\\\
x_{i+1} = x_i - \\frac {f(x_i)}{f'(x_i)},(i=0,1,2,3,...)
```

# Arguments
- `f`: `f(x)`函数
- `fd`: `f(x)`的导函数
- `x_0`: 初始估计值
- `k`: 迭代次数
- `m`: 重根树(default=1)

# Example
```jldoctest
julia> newton_raphson_method(x->x^3 +x -1, x->3x^2 +1, -0.7, 20)
0.6823278038280193
```
"""
function newton_raphson_method(f, fd, x_0, k, m=1)
    x_i = x_0
    for i = 1:k
        x_i_plus_1 = x_i - m*f(x_i)/fd(x_i)
        x_i = x_i_plus_1
    end
    return x_i
end