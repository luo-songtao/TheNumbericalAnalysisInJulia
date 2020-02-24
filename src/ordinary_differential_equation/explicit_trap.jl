"""
# 显示梯形方法(改进的欧拉方法)
    explicit_trap(f, a, b, y_0, n)
只要对欧拉公式做一个小的调整，就可以对精度有很大的提高

![](../img/explicit_trapezoid_method.png)

# Example

对比欧拉方法的结果，可以看到误差被缩小了很多

```jldoctest
julia> test_method(explicit_trap)    # 分别表示：步数 步长 误差
10×3 Array{Float64,2}:
    5.0  0.2          0.00319259 
   10.0  0.1          0.000965935
   20.0  0.05         0.000266896
   40.0  0.025        7.02111e-5 
   80.0  0.0125       1.80092e-5 
  160.0  0.00625      4.56067e-6 
  320.0  0.003125     1.14755e-6 
  640.0  0.0015625    2.87815e-7 
 1280.0  0.00078125   7.207e-8   
 2560.0  0.000390625  1.80321e-8  
```

## 显示梯形方法的局部截断误差

![](../img/explicit_trapezoid_method_local_error.png)

可见:
```math
y_{i+1} - \\omega_{i+1} = O(h^3)
```

"""
function explicit_trap(df, a, b, y_0, n)
    ω_i = y_0
    h = (b-a)/n
    for t_i in a:h:b-h    # t_i 更新后得出的是t_i+1的值，所以如果返回的是t_i处的近似值，那么应该将b-h
        g = df(t_i, ω_i)
        ω_i = ω_i + h/2 * (g + df(t_i+h, ω_i + h*g))
    end
    return ω_i
end
