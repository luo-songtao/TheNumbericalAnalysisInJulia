# 迭代方法

```@meta
CurrentModule = EquationSet
```

```@index
Pages   = ["iteration_method.md"]
```

## 雅可比方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :jacobi
```

## 高斯-赛德尔方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :gauss_seidel
```

## 连续过松弛方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :successive_over_relaxation
```
