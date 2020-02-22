# 欧拉方法

```@meta
CurrentModule = OrdinaryDifferentialEquation
```

```@index
Pages   = ["euler.md"]
```

## 欧拉方法
```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :euler
```

## Lipschitz连续

![](../img/lipschitz_continuous.png)

## 局部和全局截断误差

![](../img/local_and_global_truncation_error.png)

## 局部截断误差和全局截断误差的关系

- 局部截断误差和``h^{k+1}``成正比
- 全局截断误差和``h^k``成正比

对于欧拉方法
- 局部截断误差和``h^{2}``成正比
- 全局截断误差和``h^1``成正比

欧拉方法是一阶的

![](../img/local_and_global_truncation_error2.png)
