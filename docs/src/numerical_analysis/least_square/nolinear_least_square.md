# 非线性最小二乘

```@meta
CurrentModule = LeastSquare
```

```@index
Pages   = ["nolinear_least_square.md"]
```

## 高斯牛顿法
```@autodocs
Modules = [LeastSquare]
Filter = f -> nameof(f) == :gauss_newton
```

## Levenberg-Marquardt方法
```@autodocs
Modules = [LeastSquare]
Filter = f -> nameof(f) == :levenberg_marquardt
```
