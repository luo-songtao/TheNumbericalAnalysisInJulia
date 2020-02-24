# Runge Kutta 方法

```@meta
CurrentModule = OrdinaryDifferentialEquation
```

```@index
Pages   = ["runge_kutta.md"]
```

Runge Kutta 方法是一组ODE求解器。包含了欧拉和梯形方法，以及更复杂的高阶方法

## 中点方法

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :rk_midpoint
```

## 4阶 Runge Kutta 方法

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :rk4
```

## 计算机仿真：Hodgkin Huxley神经元
```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :hodgkin_huxley
```