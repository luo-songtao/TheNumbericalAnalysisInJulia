"""
# 霍纳法则
    horner_rule(degree, constants, x, base_points)

对多项式求值的霍纳法则(嵌套乘法)

```math
c_1 + (x-r_1)(c_2+(x-r_2)(c_3 + (x-r_3)(c_4+ (x-r_4)c_5)))
```

# Arguments
- `degree`: 多项式的阶
- `constants`: degree+1个系数数组，分别是``x^0,x^1,···,x^n``项的系数
- `x::Real`: 进行求值带入的x
- `base_points`: 每一阶的x的基点数组

# Usage
```jldoctest
julia> horner_rule(4, [-1 5 -3 3 2], 1/2, [0 0 0 0])
1.25
```
"""
function horner_rule(degree, constants, x, base_points)
    y = constants[degree+1]
    for i = degree:-1:1
        y = y.*(x-base_points[i])+constants[i]
    end
    y
end

nest = horner_rule