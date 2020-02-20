# 数值积分

```@meta
CurrentModule = NumericalCalculus
```

```@index
Pages   = ["numerical_integral.md"]
```

## 闭牛顿-科特斯方法

### (复合)梯形法则
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :trapezoidal_rule
```

### (复合)辛普森法则
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :simpson_rule
```

## 开牛顿-科特斯方法

### (复合)中点法则
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :middle_point_rule
```


## Romberg 积分
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :romberg_integral
```

## 自适应积分

### 自适应梯形法则积分
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :trap_adaptive_integral
```

### 自适应辛普森法则积分
```@autodocs
Modules = [NumericalCalculus]
Filter = f -> nameof(f) == :simpson_adaptive_integral
```
