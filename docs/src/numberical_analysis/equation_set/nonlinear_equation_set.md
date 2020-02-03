# 非线性方程组

```@meta
CurrentModule = EquationSet
```

```@index
Pages   = ["nonlinear_equation_set.md"]
```

## 多变量牛顿法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :mult_newton_method
```

## Broyden 方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :broyden
```

## BroydenⅡ 方法
```@autodocs
Modules = [EquationSet]
Filter = f -> nameof(f) == :broyden2
```