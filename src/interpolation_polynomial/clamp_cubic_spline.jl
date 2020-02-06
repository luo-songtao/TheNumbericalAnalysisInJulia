"""
# 钳制三次样条
    clamp_cubic_spline(x, y, v1, vn, k=10)

钳制三次样条和曲率调整三次样条类似，不同的是，它的一阶导数``S'_1(x_1)、S'_{n-1}(x_n)``都是可有用户选择的定义的数值。这样在样条的开始和结束端点的斜率由用户来控制。

## 约束条件
```math
S'_1(x_1) = v_1 \\\\
S'_{n-1}(x_n) = v_n
```

## 推导
因为

``\\qquad \\begin{aligned} S'_1(x_1) &= b_1 \\\\ S'_{n-1}(x_n) &= b_{n-1} + c_{n-1}\\delta_{n-1} + d_{n-1}\\delta^2_{n-1} \\end{aligned}``

可得两个新的方程：

``\\qquad \\begin{aligned} 2\\delta_1c_1 +\\delta_1c_2 &= 3(\\frac{\\Delta_1}{\\delta_1}-v_1) \\\\ \\delta_{n-1}c_{n-1} +2\\delta_{n-1}c_n &= 3(v_n-\\frac{\\Delta_{n-1}}{\\delta_{n-1}}) \\end{aligned}``
"""
function clamp_cubic_spline(x, y, v1, vn, k=10)
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

    A[1, 1:2] = [2δ[1] δ[1]]
    A[n, n-1:n] = [δ[n-1] 2δ[n-1]]
    r[1] = 3(Δ[1]/δ[1] -v1)
    r[n] = 3(vn - Δ[n-1]/δ[n-1])

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
