push!(LOAD_PATH, "../")
using OrdinaryDifferentialEquation

"""
# 向量值函数

对于常微分方程组我们需要提供的是一个向量值函数

# Example

向量值函数
对于``t\\in [0,1]``

``\\qquad \\begin{aligned} y'_1(t) &= y^2_2 - 2y_1 \\\\ y'_2 &= y-1 - y_2 - ty^2_2 \\end{aligned}``

注意每个变量都有自己的初值，如``y_1(0) = 0;y_2(0)=1``

上面的初值问题的解对应的是如下向量值函数：

``\\qquad \\begin{aligned} y_1(t) &= te^{-2t} \\\\ y_2 &= e^{-t} \\end{aligned}``

```julia
function vectorization(t, y)
    df1 = y[2]^2-2y[1]
    df2 = y[1]-y[2]-t*y[2]^2
    return Float64[df1, df2]
end
```
```jldoctest
julia> euler(vectorization, 0, 1, [0.0, 1.0], 10)    # 使用欧拉方法
2-element Array{Float64,1}:
 0.14687398022929213
 0.3643017723635318
julia> explicit_trap(vectorization, 0, 1, [0.0, 1.0], 10)    # 使用梯形方法
2-element Array{Float64,1}:
 0.13469654811260603
 0.36928227205630715
```
"""
function vectorization(t, y)
    df1 = y[2]^2-2y[1]
    df2 = y[1]-y[2]-t*y[2]^2
    return Float64[df1, df2]
end


# println(euler(ydot, 0 , 1, [0.0, 1.0], 10))
# function euler2(fs, a, b, y_0, n)
#     ω = y_0
#     h = (b-a)/n
#     for t in a:h:b-h
#         for i = 1:length(ω)
#             f = fs[i]
#             ω[i] = ω[i] + h*f(t, ω...)
#         end
#     end
#     return ω
# end