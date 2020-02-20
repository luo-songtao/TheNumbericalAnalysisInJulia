# 数值微分

```@meta
CurrentModule = NumericalCalculus
```

```@index
Pages   = ["numerical_differential.md"]
```

## 二点前向差分公式
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :two_point_forward_diff_formula
```

## 三点中心差分公式
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :three_point_mid_diff_formula
```

## 二阶导数的三点中心差分公式
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :three_point_mid_diff_formula_for_second_derivative
```

## 外推
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :extrapolation
```

