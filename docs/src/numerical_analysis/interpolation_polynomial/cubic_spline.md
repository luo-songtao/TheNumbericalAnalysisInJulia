# 三次样条插值

```@meta
CurrentModule = InterpolationPolynomial
```

```@index
Pages   = ["cubic_spline.md"]
```

## 自然三次样条
```@autodocs
Modules = [InterpolationPolynomial]
Filter = f -> nameof(f) == :natural_cubic_spline
```

## 曲率调整三次样条
```@autodocs
Modules = [InterpolationPolynomial]
Filter = f -> nameof(f) == :curvature_adjustable_cubic_spline
```

## 钳制三次样条
```@autodocs
Modules = [InterpolationPolynomial]
Filter = f -> nameof(f) == :clamp_cubic_spline
```

## 抛物线端点三次样条
```@autodocs
Modules = [InterpolationPolynomial]
Filter = f -> nameof(f) == :parabola_endpoint_cubic_spline
```

## 非纽结三次样条
```@autodocs
Modules = [InterpolationPolynomial]
Filter = f -> nameof(f) == :unknot_cubic_spline
```

## 不同三次样条方法比较

##### 下面分别是针对``x=[0 1 2 3 4 5]、y=[3 1 4 1 2 0]``进行三次样条插值的图像

```julia
x_0 = [0;1;2;3;4;5]
y_0 = [3;1;4;1;2;0]

n = length(x_0)

function generate_coordinate(coeff)
    x = []
    y = []
    for i = 1:n-1
        xs = x_0[i]:0.05:x_0[i+1]
        b, c, d = coeff[i, :]
        dx = xs.-x_0[i]
        ys = dx.*d
        ys = dx.*(ys.+c)
        ys = dx.*(ys.+b) .+ y_0[i]
        x = [x;xs]
        y = [y;ys]
    end
    return x, y
end

x1, y1 = generate_coordinate(natural_cubic_spline(x_0, y_0))
x2, y2 = generate_coordinate(curvature_adjustable_cubic_spline(x_0, y_0, 10, 10))    # 端点处曲率设置为10
x3, y3 = generate_coordinate(clamp_cubic_spline(x_0, y_0, 0, 0))    # 端点处斜率设置为0
x4, y4 = generate_coordinate(parabola_endpoint_cubic_spline(x_0, y_0))
x5, y5 = generate_coordinate(unknot_cubic_spline(x_0, y_0))

plot(
    plot!(scatter(x_0, y_0, label="base point"), x1, [y1 y2], title="natural/curvature cubic spline", label=["nature" "curvature=10"]),
    plot!(scatter(x_0, y_0, label="base point"), x3, [y3], title="clamp cubic spline", label="clamp"),
    plot!(scatter(x_0, y_0, label="base point"), x4, [y4], title="parabola endpoint cubic spline", label="parabola"),
    plot!(scatter(x_0, y_0, label="base point"), x5, [y5], title="unknot cubic spline", label="unknot")
)
```

![](../img/cubic_spline.png)
