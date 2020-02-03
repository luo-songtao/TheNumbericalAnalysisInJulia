module EquationSet

include("equation_set/naive_gauss_elimination.jl")
export naive_gauss_elimination

include("equation_set/lu_factorization.jl")
export lu_factorization

include("equation_set/lu_factorization_with_permutation.jl")
export lu_factorization_with_permutation

include("equation_set/jacobi.jl")
export jacobi

include("equation_set/gauss_seidel.jl")
export gauss_seidel

include("equation_set/sor.jl")
export sor

include("equation_set/cholesky_decomposition.jl")
export cholesky_decomposition

include("equation_set/conjugate_gradient.jl")
export conjugate_gradient

include("equation_set/conjugate_gradient_with_pre_condition.jl")
export conjugate_gradient_with_pre_condition

include("equation_set/ssor.jl")
export ssor

include("equation_set/multiple_newton_method.jl")
export mult_newton_method

include("equation_set/broyden.jl")
export broyden

include("equation_set/broyden2.jl")
export broyden2

end