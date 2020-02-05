# 多项式求值

```@meta
CurrentModule = Polyomial
```

## Content
```@contents
Pages = map(file -> joinpath("equation_solving", file), readdir("equation_solving"))
Depth = 3
```

## Index
```@index
Pages   = ["polyomial.md"]
```

#### 霍纳方法
```@autodocs
Modules = [Polyomial]
Filter = f -> nameof(f) == :horner_rule
```