# 用于对称正定矩阵的方法

```@meta
CurrentModule = EquationSet
```

```@index
Pages   = ["methods_about_SPD_matrix.md"]
```

## 楚列斯基分解法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :cholesky_decomposition 
```

## 共轭梯度方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :conjugate_gradient 
```

## 预条件共轭梯度方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :conjugate_gradient_with_pre_condition 
```
