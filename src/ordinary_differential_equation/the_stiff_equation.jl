push!(LOAD_PATH, "../")
using OrdinaryDifferentialEquation

df(t, y) = 10*(y.-1)

# rk45_dp(df, 0, 100, 0.01, [0.5], 1, 1e-4, 1e-4, 0.9, 1.5, 0.5)

"""
对于刚性方程，使用如Runge Kutta: Dormand Prince 4/5这样的方法，很难收敛到正确的解，而且这类方法都是显示的。

对于刚性方程，
"""
function example_backward_euler(a, b, n, y_0)
    ω_i = y_0
    h = (b-a)/n
    for t_i in a:h:b-h    # t_i 更新后得出的是t_i+1的值，所以如果返回的是t_i处的近似值，那么应该将b-h
        ω_i = (ω_i + 10h)/(1+10h)
        println(ω_i)
    end
    return ω_i
end

# example_backward_euler(0, 100, 100, 0.5)