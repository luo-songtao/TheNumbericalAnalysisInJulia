"""
# 二分法

    bisect(fn, a, b, err)
二分法求解方程解

# Arguments
- `fn`: 函数
- `a`: 左区间点
- `b`: 右区间点
- `err`: 误差

# Example
```jldoctest
julia> bisect(x->x^3+x-1, 0, 1, 0.00005)
0.682342529296875
```
"""
function bisect(fn, a ,b, err)
    if sign(fn(a)) * sign(fn(b)) >= 0
        error("f(a)f(b) < 0 not satisfied!")
    end
    fa = fn(a)
    fb = fn(b)

    while (b-a)/2 > err
        c = (a+b)/2
        fc = fn(c)
        if fc == 0
            return c
        elseif sign(fc)*sign(fa) < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
    end
    return (a+b)/2    # 返回新的中点作为最优估计
end