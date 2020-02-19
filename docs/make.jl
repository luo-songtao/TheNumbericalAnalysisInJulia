push!(LOAD_PATH,"../src/")
using Documenter, EquationSolving, EquationSet, InterpolationPolynomial, LeastSquare

modules = [EquationSolving, EquationSet, InterpolationPolynomial, LeastSquare]

DocMeta.setdocmeta!(EquationSolving, :DocTestSetup, :(using EquationSolving); recursive=true)
DocMeta.setdocmeta!(EquationSet, :DocTestSetup, :(using EquationSet); recursive=true)
DocMeta.setdocmeta!(InterpolationPolynomial, :DocTestSetup, :(using InterpolationPolynomial;using EquationSet); recursive=true)
DocMeta.setdocmeta!(LeastSquare, :DocTestSetup, :(using LeastSquare;); recursive=true)

makedocs(
    modules = modules,
    sitename = "数值分析(Julia语言描述)",
    doctest = true,
    format = Documenter.HTML(assets = [
        asset("assets/logo.png", class=:ico, islocal=true)
    ]),
    pages = [
        "Home"=>"index.md",
        "求解一元方程" => [
            "目录" => "numerical_analysis/equation_solving.md",
            "numerical_analysis/equation_solving/bisect.md",
            "numerical_analysis/equation_solving/fixed_point_iteration.md",
            "numerical_analysis/equation_solving/newton_method.md",
            "numerical_analysis/equation_solving/secant.md"
        ],
        "求解n元方程组" => [
            "目录" => "numerical_analysis/equation_set.md",
            "numerical_analysis/equation_set/gauss_elimination.md",
            "numerical_analysis/equation_set/lu_factorization.md",
            "numerical_analysis/equation_set/iteration_method.md",
            "numerical_analysis/equation_set/methods_about_SPD_matrix.md",
            "numerical_analysis/equation_set/nonlinear_equation_set.md"
        ],
        "插值多项式" => [
            "目录" => "numerical_analysis/interpolation_polynomial.md",
            "numerical_analysis/interpolation_polynomial/newton_difference_quotient.md",
            "numerical_analysis/interpolation_polynomial/error_and_runge_phenomenon.md",
            "numerical_analysis/interpolation_polynomial/chebyshev_interpolation.md", 
            "numerical_analysis/interpolation_polynomial/cubic_spline.md",
            "numerical_analysis/interpolation_polynomial/bezier_curve.md",
        ],
        "最小二乘" => [
            "目录" => "numerical_analysis/least_square.md",
            "numerical_analysis/least_square/normal_equation.md",
            "numerical_analysis/least_square/gram_schmidt_orthogon.md",
            "numerical_analysis/least_square/least_square_and_qr.md",
            "numerical_analysis/least_square/householder_reflector.md",
            "numerical_analysis/least_square/gmres.md",
            "numerical_analysis/least_square/nolinear_least_square.md",
        ]
    ]
)

deploydocs(
    repo = "github.com/luo-songtao/TheNumericalAnalysisInJulia.git",
    devbranch = "dev",
    devurl = "dev"
)