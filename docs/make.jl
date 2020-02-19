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
            "目录" => "numberical_analysis/equation_solving.md",
            "numberical_analysis/equation_solving/bisect.md",
            "numberical_analysis/equation_solving/fixed_point_iteration.md",
            "numberical_analysis/equation_solving/newton_method.md",
            "numberical_analysis/equation_solving/secant.md"
        ],
        "求解n元方程组" => [
            "目录" => "numberical_analysis/equation_set.md",
            "numberical_analysis/equation_set/gauss_elimination.md",
            "numberical_analysis/equation_set/lu_factorization.md",
            "numberical_analysis/equation_set/iteration_method.md",
            "numberical_analysis/equation_set/methods_about_SPD_matrix.md",
            "numberical_analysis/equation_set/nonlinear_equation_set.md"
        ],
        "插值多项式" => [
            "目录" => "numberical_analysis/interpolation_polynomial.md",
            "numberical_analysis/interpolation_polynomial/newton_difference_quotient.md",
            "numberical_analysis/interpolation_polynomial/error_and_runge_phenomenon.md",
            "numberical_analysis/interpolation_polynomial/chebyshev_interpolation.md", 
            "numberical_analysis/interpolation_polynomial/cubic_spline.md",
            "numberical_analysis/interpolation_polynomial/bezier_curve.md",
        ],
        "最小二乘" => [
            "目录" => "numberical_analysis/least_square.md",
            "numberical_analysis/least_square/normal_equation.md",
            "numberical_analysis/least_square/gram_schmidt_orthogon.md",
            "numberical_analysis/least_square/least_square_and_qr.md",
            "numberical_analysis/least_square/householder_reflector.md",
            "numberical_analysis/least_square/gmres.md",
            "numberical_analysis/least_square/nolinear_least_square.md",
        ]
    ]
)

deploydocs(
    repo = "github.com/luo-songtao/TheNumbericalAnalysisInJulia.git",
    devbranch = "dev",
    devurl = "dev"
)