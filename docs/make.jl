push!(LOAD_PATH,"../src/")
using Documenter, EquationSolving, EquationSet, Polyomial

modules = [EquationSolving, EquationSet, Polyomial]

DocMeta.setdocmeta!(EquationSolving, :DocTestSetup, :(using EquationSolving); recursive=true)
DocMeta.setdocmeta!(EquationSet, :DocTestSetup, :(using EquationSet); recursive=true)
DocMeta.setdocmeta!(Polyomial, :DocTestSetup, :(using Polyomial); recursive=true)

makedocs(
    modules = modules,
    sitename = "数值分析--Julia实现",
    doctest = true,
    pages = [
        "目录"=>"index.md",
        "求解方程"=> "numberical_analysis/equation_solving.md",
        "方程组"=> "numberical_analysis/equation_set.md",
        "多项式求值" => "numberical_analysis/polyomial.md"
    ]
)
