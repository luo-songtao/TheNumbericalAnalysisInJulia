
module Polyomial

export horner_rule

"""
# 霍纳法则
    horner_rule(rank, constants, x, base_points)

对多项式求值的霍纳法则(嵌套乘法)
```math
c_1 + (x-r_1)(c_2+(x-r_2)(c_3 + (x-r_3)(c_4+ (x-r_4)c_5)))
```

# Arguments
- `rank`: 多项式的阶
- `constants`: rank+1个系数数组，第一个是常数项
- `x::Real`: 进行求值带入的x
- `base_points`: 每一阶的x的基点数组

# Usage
```jldoctest
julia> horner_rule(4, [-1 5 -3 3 2], 1/2, [0 0 0 0])
1.25
```
"""
function horner_rule(rank, constants, x, base_points)
    y = constants[rank+1]
    for i = rank:-1:1
        y = y.*(x-base_points[i])+constants[i]
    end
    y
end

end
