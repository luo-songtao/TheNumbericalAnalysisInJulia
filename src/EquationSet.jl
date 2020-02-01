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

include("equation_set/successive_over_relaxation.jl")
export successive_over_relaxation

end