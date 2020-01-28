module EquationSolving

include("equation_solving/bisect.jl")
export bisect

include("equation_solving/fixed_point_iteration.jl")
export fixed_point_iteration

include("equation_solving/newton_raphson_method.jl")
export newton_raphson_method

end