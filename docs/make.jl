push!(LOAD_PATH,"../src/")
using Documenter, EquationSolving, Polyomial

modules = [EquationSolving, Polyomial]

DocMeta.setdocmeta!(EquationSolving, :DocTestSetup, :(using EquationSolving); recursive=true)
DocMeta.setdocmeta!(Polyomial, :DocTestSetup, :(using Polyomial); recursive=true)

makedocs(
    modules = modules,
    sitename = "数值分析--Julia实现",
    doctest = true,
    pages = [
        "目录"=>"index.md",
        "求解方程"=> "numberical_analysis/equation_solving.md",
        "多项式求值" => "numberical_analysis/polyomial.md"
    ]
)
