"""
# 非纽结三次样条
    unknot_cubic_spline(x, y, k=10)

## 约束条件

``d_1 = d_2, d_{n-2}=d_{n-1}``

等价于: ``S'''_1{x_2} = S'''_2(x_2), S'''_{n-2}(x_{n-1})=S'''_{n-1}(x_{n-1})``

可以看到``S_1``和``S_2``在``x_2``的0、1、2、3阶导数均相等，意味着不再将``x_2``作为基点，在``[x_1, x_3]``上``S_1=S_2``;同理``x_{n-1}``不再作为基点，``S_{n-2}=S_{n-1}``

## 推导

``\\qquad \\begin{aligned} d_1=d_2 &\\Rightarrow \\frac{c_2-c_1}{\\delta_1} = \\frac{c_3-c_2}{\\delta_2} \\\\  &\\Rightarrow \\delta_2c_1 - (\\delta_1+\\delta_2)c_2 + \\delta_1c_3=0 \\\\ d_{n-2} = d_{n-1} &\\Rightarrow \\delta_{n-1}c_{n-2} - (\\delta_{n-2} + \\delta_{n-1})c_{n-1} + \\delta_{n-2}c_n = 0 \\end{aligned}``
"""
function unknot_cubic_spline(x, y, k=10)
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

    A[1, 1:3] = [δ[2] -(δ[1]+δ[2]) δ[1]]
    A[n, n-2:n] = [δ[n-1] -(δ[n-2]+δ[n-1]) δ[n-2]]

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


