"""
# Levenberg Marquardt 方法
    levenberg_marquardt(r, Dr, x_0, k, λ)
对于非线性最小二乘问题中，如果定义好的模型在计算得到了条件数比较差的雅可比矩阵的话，往往求出的解会有比较大的误差，因此
Levenberg Marquardt 方法使用了**正则化**来修复这个问题

其中用``\\lambda diag J^T_rJ_r``来强化对角线元素的作用，以改善条件数

```math
(J^T_rJ_r + \\lambda diag J^T_rJ_r) \\Delta{x_k} = -J^T_rr(x_k)
```

当``\\lambda = 0``时，该方法就是高斯牛顿方法。另外``\\lambda``通常看作一个常数，但该方法中常常使用不同的``\\lambda``以适应问题，一般策略是：
- 只要余下的平方误差和在每步降低，那么就使用一个因子降低``\\lambda`` 
- 如果误差升高，则反之，使用因子升高``\\lambda``

# Example

使用Levenberg Marquardt将模型``y = c_1 e^{-c_2(x-c_3)^2}``拟合到数据点(1,3),(2,5),(2,7),(3,5),(4,1)

```julia
function test_levenberg_marquardt()
    points = [
        (1,3)
        (2,5)
        (2,7)
        (3,5)
        (4,1)
    ]

    r = [(c1,c2,c3) -> c1*exp(1)^(-c2*(x-c3)^2)-y for (x,y) in points]
    f_c1 = [(c1,c2,c3)-> exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    f_c2 = [(c1,c2,c3)-> -c1*(x-c3)^2*exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    f_c3 = [(c1,c2,c3)-> 2c1*c2*(x-c3)*exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    Dr = hcat([f_c1,f_c2,f_c3]...)
    x = round.(levenberg_marquardt(r, Dr, [1;1;1], 1000, 50), digits=6)
    # x = round.(levenberg_marquardt(r, Dr, [1;1;1], 1000, 0), digits=6)   无法收敛
    @assert reshape(x, 1,3) == [6.300046 0.508648 2.248735]
    return x
end
```
λ固定50，迭代1000步的收敛结果：
```jldoctest
julia> test_levenberg_marquardt()
3×1 Array{Float64,2}:
 6.300046
 0.508648
 2.248735
```
"""
function levenberg_marquardt(r, Dr, x_0, k, λ)
    m,n = size(Dr)
    x_k = x_0
    for k = 1:k
        r_k = Float64[r_f(x_k...) for r_f in r]
        # println("RMSE: ", sqrt(sum(r_k.^2)))
        J_k = Float64[dr_f(x_k...) for dr_f in Dr]
        H = J_k'J_k
        diag_H = diagm(diag(H))
        Δx = naive_gauss_elimination((H+λ*diag_H), -J_k'r_k, n)
        x_k = x_k + Δx
    end
    return x_k
end

push!(LOAD_PATH,"../")
using LinearAlgebra
using EquationSet

function test_levenberg_marquardt()
    points = [
        (1,3)
        (2,5)
        (2,7)
        (3,5)
        (4,1)
    ]

    r = [(c1,c2,c3) -> c1*exp(1)^(-c2*(x-c3)^2)-y for (x,y) in points]
    f_c1 = [(c1,c2,c3)-> exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    f_c2 = [(c1,c2,c3)-> -c1*(x-c3)^2*exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    f_c3 = [(c1,c2,c3)-> 2c1*c2*(x-c3)*exp(1)^(-c2*(x-c3)^2) for (x,y) in points]
    Dr = hcat([f_c1,f_c2,f_c3]...)
    x = round.(levenberg_marquardt(r, Dr, [1;1;1], 1000, 50), digits=6)
    # x = round.(levenberg_marquardt(r, Dr, [1;1;1], 1000, 0), digits=6)   无法收敛
    @assert reshape(x, 1,3) == [6.300046 0.508648 2.248735]
    return x
end
