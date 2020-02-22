# 常微分方程组

```@meta
CurrentModule = OrdinaryDifferentialEquation
```

```@index
Pages   = ["equations.md"]
```

对于微分方程组的近似求解，可以通过单个微分方程近似方法的简单扩展得到。

微分方程的阶指的是出现在方程中的最高阶导数。

## 向量值函数

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :vectorization
```

## 高阶方程处理方式

单个高阶微分方程可以转化为一个方程组

![](../img/high_order_equation.png)