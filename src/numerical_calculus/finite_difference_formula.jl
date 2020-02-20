"""
# 二点前向差分公式
    two_point_forward_diff_formula(f, x, h)
二点前向差分公式是近似一阶导数的二阶近似，误差为O(h)

由定义f在x点的导数
```math
f(x) = \\lim_{h\\Rightarrow 0} \\frac {f(x+h)-f(x)}{h}
```

如果f是二阶连续可微函数，则:
```math
f(x+h) = f(x) + hf'(x) + \\frac 12 h^2f''(c)
```

其中c在x和x+h之间

则有:
```math
f'(x) = \\frac {f(x+h)-f(x)}{h} - \\frac 12 hf''(c)
```

# Example
```jldoctest
julia> a = two_point_forward_diff_formula(x->1/x, 2, 0.1)
-0.23809523809523836
julia> e = a - (-2^(-2))
0.01190476190476164
```
"""
function two_point_forward_diff_formula(f, x, h)
    return (f(x+h) - f(x))/h
end


"""
# 三点中心差分公式
    three_point_mid_diff_formula(f, x, h)
三点中心差分公式是近似一阶导数的二阶近似，误差为``O(h^2)``，

由定义f在x点的导数
```math
f(x) = \\lim_{h\\Rightarrow 0} \\frac {f(x+h)-f(x)}{h}
```

如果f是三阶连续可微函数，则:
```math
f(x+h) = f(x) + hf'(x) + \\frac 12 h^2f''(x) + \\frac 16 h^3f'''(c_1) \\\\
f(x-h) = f(x) - hf'(x) + \\frac 12 h^2f''(x) - \\frac 16 h^3f'''(c_2)
```

其中c在x-h和x+h之间

两式相减, 并取相同的c，则有:
```math
f'(x) = \\frac {f(x+h)-f(x-h)}{2h} - \\frac {1}{6} h^2f'''(c)
```

# Example
```jldoctest
julia> a = three_point_mid_diff_formula(x->1/x, 2, 0.1)
-0.2506265664160401
julia> e = a - (-2^(-2))
-0.0006265664160400863
```
"""
function three_point_mid_diff_formula(f, x, h)
    return (f(x+h) - f(x-h))/(2h)
end

"""
# 二阶导数的三点中心差分公式
    three_point_mid_diff_formula_for_second_derivative(f, x, h)
三点中心差分公式是近似二阶导数的二阶近似，误差为``O(h^2)``

由定义f在x点的导数
```math
f(x) = \\lim_{h\\Rightarrow 0} \\frac {f(x+h)-f(x)}{h}
```

如果f是三阶连续可微函数，则:
```math
f(x+h) = f(x) + hf'(x) + \\frac 12 h^2f''(x) + \\frac 16 h^3f'''(c_1) \\\\
f(x-h) = f(x) - hf'(x) + \\frac 12 h^2f''(x) - \\frac 16 h^3f'''(c_2)
```

其中c在x-h和x+h之间

两式加, 并取相同的c，则有:
```math
f'(x) = \\frac {f(x+h)-2f(x)+f(x-h)}{h^2} - \\frac {1}{12} h^2f^{(4)}(c)
```

# Example
```jldoctest
julia> a = three_point_mid_diff_formula_for_second_derivative(x->1/x, 2, 0.1)
0.2506265664160345
julia> e = a - (2*2^(-3))
0.0006265664160344797
```
"""
function three_point_mid_diff_formula_for_second_derivative(f, x, h)
    return (f(x+h) - 2f(x) + f(x-h))/(h^2)
end

"""
# 外推公式

假设有n阶公式``F(h)``，近似一个给定量Q，这个阶数意味着：
```math
Q \\approx F(h) + Kh^n
```
K大约是h区间上的一个常数

但如果在公式中使用``\\frac h2``替代公式中的h，那么误差就由常数K乘上``h^n``变为常数乘上``(\\frac h2)^2``，或者说以因子``2^n``在降低，也就是说，我们期望：
```math
Q - F(\\frac h2) \\approx \\frac {Q-F(h)}{2^n}
```

这样得到一个``F(h)``的**n阶外推公式**，也称为理查德森外推：

```math
Q \\approx \\frac {2^nF(\\frac h2)- F(h)}{2^n-1}
```
"""
function extrapolation()

end