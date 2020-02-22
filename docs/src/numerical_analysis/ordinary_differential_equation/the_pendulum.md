# 计算机仿真: 钟摆模型

```@meta
CurrentModule = OrdinaryDifferentialEquation
```

```@index
Pages   = ["the_pendulum.md"]
```

## 钟摆模型

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :the_pendulum
```

## 无衰减钟摆模型

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :the_undamped_pendulum_model
```

## 衰减钟摆模型

```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :the_damped_pendulum_model
```

## 受力衰减钟摆模型
```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :the_forced_damped_pendulum_model
```

## 双钟摆模型
```@autodocs
Modules = [OrdinaryDifferentialEquation]
Filter = f -> nameof(f) == :the_double_pendulum_model
```
