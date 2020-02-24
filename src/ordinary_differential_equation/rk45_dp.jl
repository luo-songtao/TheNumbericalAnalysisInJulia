using LinearAlgebra

# Dormand-Prince coefficients
const rk_dp_coefficients = Rational{Int64}[
        0        0            0                   0            0      0 0
    1//5        0            0                   0            0      0 0
    3//40       9//40        0                   0            0      0 0
    44//45      -56//15      32//9               0            0      0 0
    19372//6561 -25360//2187 64448//6561 -212//729            0      0 0
    9017//3168  -355//33     46732//5247   49//176 -5103//18656      0 0
    35//384     0            500//1113    125//192  -2187//6784 11//84 0],
    [35//384         0     500//1113       125//192      -2187//6784         11//84      0
    5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
    [0, 1//5, 3//10, 4//5, 8//9, 1, 1]

"""
# Runge Kutta: Dormand Prince 4/5
    rk45_dp(df, t_start, t_end, h_start, y_0, n, abs_tol=1e-6, rel_tol=1e-6, fac=0.8, facmax=2, facmin=0.6)
Runge Kutta: Dormand Prince 4阶/5阶嵌入对

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

## Dormand Prince 系数公式

![](../img/rk_dp45.png)

## 步长控制

![](../img/automatic_step_size_control.png)

# Example

```julia
function test_rk45_dp()
    df(t, y) = y.*t .+ t^3
    rk45_dp(df, 0, 1, 0.01, [1], 1)

    println()
    
    function dfs(t, y)
        df1 = y[2]^2-2y[1]
        df2 = y[1]-y[2]-t*y[2]^2
        return Float64[df1, df2]
    end
    rk45_dp(dfs, 0, 1, 0.01, [0 1], 2)
end
```

```julia
julia> test_rk45_dp()
1: 0.01  [1.0000500037500626]
2: 0.03  [1.0004503037955677]
3: 0.07  [1.002459011107573]
4: 0.15000000000000002  [1.0114405576720515]
5: 0.31000000000000005  [1.0515693458693955]
6: 0.626908238067097  [1.258417618751055]
7: 0.8831300641373249  [1.650843872358512]
8: 1.0  [1.946164140264331]

1: 0.01  [0.00980198673299806 0.9900498337492146]
2: 0.03  [0.028252936002962036 0.9704455335515371]
3: 0.07  [0.060855076181021095 0.9323938201001817]
4: 0.15000000000000002  [0.11112271393407858 0.8607079888957403]
5: 0.24842403330400986  [0.15115239787448706 0.7800291612233432]
6: 0.35023793315011253  [0.17384016930556773 0.7045205113148135]
7: 0.4570032978500713  [0.18321922428238568 0.6331783401091992]
8: 0.569316454560961  [0.18232719398790645 0.5659122339459881]
9: 0.6879863459066563  [0.1737809427101828 0.5025871963036914]
10: 0.8140215103284122  [0.1598031853267181 0.4430727708639373]
11: 0.9487341496417322  [0.1422604382262948 0.387231004857348]
12: 1.0  [0.1353351750388877 0.36787954090370917]]
```

"""
function rk45_dp(df, t_start, t_end, h_start, y_0, n, abs_tol=1e-6, rel_tol=1e-6, fac=0.8, facmax=2, facmin=0.6)
    coeff = rk_dp_coefficients
    t_i = t_start
    h_i = h_start
    w_i = y_0
    step = 0
    while t_i < t_end
        if t_i+h_i > t_end
            h_i = t_end - t_i
        end
        s = zeros(7,n)
        for k = 1:6
            s[k,:] = df(t_i+coeff[3][k]*h_i, w_i + h_i*coeff[1][k,:]'*s)
        end
        z_ip1 = w_i+h_i*coeff[2][1,:]'*s
        s[7, :] = df(t_i+h_i, z_ip1)
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

function test_rk45_dp()
    df(t, y) = y.*t .+ t^3
    rk45_dp(df, 0, 1, 0.01, [1], 1)
    println()
    function dfs(t, y)
        df1 = y[2]^2-2y[1]
        df2 = y[1]-y[2]-t*y[2]^2
        return Float64[df1, df2]
    end
    rk45_dp(dfs, 0, 1, 0.01, [0 1], 2)
end

# test_rk45_dp()