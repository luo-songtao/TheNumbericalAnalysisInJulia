"""
# 曲率调整三次样条
    curvature_adjustable_cubic_spline(x, y, r1, rn, k=10)
曲率调整三次样条和自然三次样条的不同之处在于，其中的``S''_1(x_1)、S''_{n-1}(x_n)``都是可由用户选择的任意数值，而不在是固定的0。这样在样条的开始和结束端点的曲率由用户来控制。

```math
S''_1(x_1) = 2c_1 = r_1 \\\\
S''_{n-1}(x_n) = 2c_n = r_n
```
"""
function curvature_adjustable_cubic_spline(x, y, r1, rn, k=10)
    n = length(x)
    # Ac = r
    A = zeros(n,n)
    c_0 = zeros(n)
    r = zeros(n)

    δ = zeros(n-1)
    Δ = zeros(n-1)

    for i = 1:n-1
        δ[i] = x[i+1] - x[i]
        Δ[i] = y[i+1] - y[i]
    end

    for i = 2:n-1
        A[i, i] = 2(δ[i-1] + δ[i])
        A[i, i-1] = δ[i-1]
        A[i, i+1] = δ[i]
        r[i] = 3(Δ[i]/δ[i] - Δ[i-1]/δ[i-1])
    end

    A[1, 1] = 2
    A[n, n] = 2
    r[1] = r1
    r[n] = rn
    c = gauss_seidel(A, r, c_0, n, k)    # 使用高斯-塞德尔方法求解c

    b = zeros(n)
    d = zeros(n)
    for i = 1:n-1
        b[i] = (c[i+1] - c[i]) / (3δ[i])
        d[i] = Δ[i]/δ[i] - (δ[i]/3)*(2c[i] + c[i+1])
    end
    return [b c d]
end

push!(LOAD_PATH,"../")
using EquationSet
