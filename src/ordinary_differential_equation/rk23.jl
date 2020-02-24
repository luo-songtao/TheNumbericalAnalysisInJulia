using LinearAlgebra

# Dormand-Prince coefficients
const rk23_coefficients = Rational{Int64}[
                                            0        0        0 
                                            1        0        0 
                                         1//4     1//4        0],
                                        [1//6     1//6     2//3
                                         1//2     1//2        0],
                                        [0, 1, 1//2]

"""
# Runge Kutta 2/3
    rk23(df, t_start, t_end, h_start, y_0, n, abs_tol=1e-6, rel_tol=1e-6, fac=0.8, facmax=2, facmin=0.6)
Runge Kutta 2阶/3阶嵌入对

# Arguments
- `df`: 常微分方程组函数
- `t_start`: 区间左端点
- `t_end`: 区间右端点
- `h_start`: 初始步长
- `y_0`: 方程组的初值条件
- `n`: 方程组数量
- `abs_tol`: 绝对误差（default=1e-6）
- `rel_tol`: 相对误差（default=1e-6）
- `fac`: 安全因子，控制步长的变化程度（default=0.8）
- `facmax`: 最大安全因子（default=2）
- `facmin`: 最小安全因子（default=0.6）

## Runge Kutta 2阶/3阶嵌入对系数公式

![](../img/rk2_3.png)

## 步长控制

同`rk45_dp`

# Example

```julia
function test_rk23()
    df(t, y) = y.*t .+ t^3
    rk23(df, 0, 1, 0.01, [1], 1)

    println()

    function dfs(t, y)
        df1 = y[2]^2-2y[1]
        df2 = y[1]-y[2]-t*y[2]^2
        return Float64[df1, df2]
    end
    rk23(dfs, 0, 1, 0.01, [0 1], 2)
end
```

```julia
julia> test_rk23()
1: 0.01  [1.0000500033334168]
2: 0.03  [1.000450296719085]
3: 0.061595012215853465  [1.0019023255625061]
4: 0.0933735421521446  [1.0043877656983238]
......
49: 0.9524037994111383  [1.8145313339944313]
50: 0.9711874755179333  [1.864475164282171]
51: 0.9882896060971613  [1.9121982307523833]
52: 1.0  [1.946163064155229]

1: 0.01  [0.009802001604500815 0.9900498274919775]
2: 0.01808072705229334  [0.0174386027842691 0.9820817398768917]
3: 0.025533702395294917  [0.024262521910342574 0.9747895152927601]
4: 0.03277497389344002  [0.03069550399294721 0.9677562937426204]
......
98: 0.966926093345866  [0.13980806669852433 0.3802500259268951]
99: 0.981516365338917  [0.13783628070487436 0.374742350372051]
100: 0.9962540822738747  [0.1358423218572797 0.36926000038324647]
101: 1.0  [0.13533537093325995 0.3678793707911274]
```

"""
function rk23(df, t_start, t_end, h_start, y_0, n, abs_tol=1e-6, rel_tol=1e-6, fac=0.8, facmax=2, facmin=0.6)
    coeff = rk23_coefficients
    t_i = t_start
    h_i = h_start
    w_i = y_0
    step = 0
    while t_i < t_end
        if t_i+h_i > t_end
            h_i = t_end - t_i
        end
        s = zeros(3,n)
        for k = 1:3
            s[k,:] = df(t_i+coeff[3][k]*h_i, w_i + h_i*coeff[1][k,:]'*s)
        end
        z_ip1 = w_i+h_i*coeff[2][1,:]'*s
        w_ip1 = w_i+h_i*coeff[2][2,:]'*s

        sc = abs_tol + max(norm(z_ip1), norm(w_i))*rel_tol
        err = norm((z_ip1 - w_ip1)./sc)
        new_h_i = h_i*min(facmax, max(facmin, fac*(1/err)^(1/5) ) )
        if err > 1
            continue
        end
        w_i = z_ip1
        t_i = t_i + h_i
        h_i = new_h_i
        step = step + 1
        println("$step: $t_i  $w_i")
    end
end

function test_rk23()
    df(t, y) = y.*t .+ t^3
    rk23(df, 0, 1, 0.01, [1], 1)
    println()
    function dfs(t, y)
        df1 = y[2]^2-2y[1]
        df2 = y[1]-y[2]-t*y[2]^2
        return Float64[df1, df2]
    end
    rk23(dfs, 0, 1, 0.01, [0 1], 2)
end

# test_rk23()