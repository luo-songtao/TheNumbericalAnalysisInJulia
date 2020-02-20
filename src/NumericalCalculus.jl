module NumericalCalculus

include("numerical_calculus/finite_difference_formula.jl")
export two_point_forward_diff_formula, three_point_mid_diff_formula, three_point_mid_diff_formula_for_second_derivative
export extrapolation

include("numerical_calculus/trapezoidal_rule.jl")
export trapezoidal_rule

include("numerical_calculus/simpson_rule.jl")
export simpson_rule

include("numerical_calculus/middle_point_rule.jl")
export middle_point_rule

include("numerical_calculus/romberg_integral.jl")
export romberg_integral

include("numerical_calculus/adaptive_integral.jl")
export trap_adaptive_integral, simpson_adaptive_integral

end