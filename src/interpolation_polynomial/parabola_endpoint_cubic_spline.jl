"""
# 抛物线端点三次样条
    parabola_endpoint_cubic_spline(x, y, k=10)

## 约束条件

``d_1=0=d_{n-1}``

通过定义``d_1=0=d_{n-1}``使得样条的起始和结束部分的``S_1``和``S_{n-1}``至多2阶。

可通过要求``c_1=c_2,c_{n-1}=c_n``使得约束条件成立
"""
function parabola_endpoint_cubic_spline(x, y, k=10)
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

    A[1, 1:2] = [1 -1]
    A[n, n-1:n] = [1 -1]

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
