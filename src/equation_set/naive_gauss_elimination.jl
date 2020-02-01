"""
# 朴素的高斯消去法
    naive_gauss_elimination(a, b, n)
朴素的高斯消去法求解方程组。高斯方法是求解线性方程组的直接方法。

1. 消去步骤
    - 注意：消去过程采取j+1:n的方式迭代，完成后应将a[i,j]置为0，
      但这里并没有这么做，是因为在后面的步骤中，并没有用到它，所以省去了置为0的过程
2. 回代步骤
    - 从底部开始，逐渐向上求解对应的方程

注意：当前实现的算法如果遇到0主元会抛出异常并终止

```math
Ax=b
```

# Arguments
- `a`: 表示系数矩阵A
- `b`: 表示常数项b
- `n`: 方程数

# Usage
```jldoctest
julia> m = [1 2 -1;2 1 -2;-3 1 1]
3×3 Array{Int64,2}:
  1  2  -1
  2  1  -2
 -3  1   1
julia> b = [3; 3; -6]
3-element Array{Int64,1}:
  3
  3
 -6
julia> x = naive_gauss_elimination(m,b,3)
3×1 Array{Float64,2}:
 3.0
 1.0
 2.0
julia> m    # 其实表示的应是[1 2 -1; 0 -3 0; 0 0 -2]
3×3 Array{Int64,2}:
  1   2  -1
  2  -3   0
 -3   7  -2
julia> b
3-element Array{Int64,1}:
  3
 -3
 -4
```
"""
function naive_gauss_elimination(a, b, n)
    # 1. 消去步骤
    # 注意：消去过程采取j+1:n的方式迭代，完成后应将a[i,j]置为0，
    # 但这里并没有这么做，是因为在后面的步骤中，并没有用到它，所以省去了置为0的过程
    for j = 1:n-1
        if abs(a[j,j]) < eps(1.0)    # 等价于2.0^-52
            error("Zero pivot encounted")
        end
        for i = j+1:n
            mult = a[i,j]/a[j,j]    # 行变换乘子
            for k = j+1:n
                a[i,k] = a[i,k] - mult * a[j,k]
            end
            b[i] = b[i] - mult*b[j]
        end
    end
    x = zeros(n, 1)
    # 2. 回代步骤
    for i = n:-1:1
        for j = i+1:n
            b[i] = b[i] - a[i,j]*x[j]
        end
        x[i] = b[i]/a[i,i]
    end
    return x
end
