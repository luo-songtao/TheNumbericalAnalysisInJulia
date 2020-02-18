# 最小二乘与QR分解

```@meta
CurrentModule = LeastSquare
```

```@index
Pages   = ["least_square_and_qr.md"]
```

## 不完全的QR分解
```@autodocs
Modules = [LeastSquare]
Filter = f -> nameof(f) == :incomplete_qr
```

## 完全QR分解
```@autodocs
Modules = [LeastSquare]
Filter = f -> nameof(f) == :complete_qr
```

## 通过QR分解实现最小二乘
```@autodocs
Modules = [LeastSquare]
Filter = f -> nameof(f) == :least_square_by_complete_qr
```
