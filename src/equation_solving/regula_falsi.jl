"""
# 试位方法
    regula_falsi(f, a, b, k)
试位方法是割线方法推广，类似二分法，但其中的中点被类似割线方法的近似所替换

## 试位方法迭代公式
``\\begin{aligned} a, b &= 初始估计值\\\\ c &= a - \\frac {f(a)(a-b)}{f(a)-f(b)} = \\frac {bf(a)-af(b)}{f(a)-f(b)} \\end{aligned}``

# Arguments
- `f`: `f(x)`函数
- `a`: 初始估计值1
- `b`: 初始估计值2
- `k`: 迭代次数

# Example
```jldoctest
julia> regula_falsi(x->x^3 +x -1, 0, 1, 100)
0.6823278038280193
```
"""
function regula_falsi(f, a, b, k)
    if sign(f(a)) * sign(f(b)) >= 0
        error("f(a)f(b) < 0 not satisfied!")
    end
    fa = f(a)
    fb = f(b)
    
    for i = 1:k
        c = (b*fa-a*fb)/(fa-fb)
        fc = f(c)
        if fc == 0; return c; end
        if fa*fc < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
        if fa == fb; return a; end
    end

    return (b*fa-a*fb)/(fa-fb)
end