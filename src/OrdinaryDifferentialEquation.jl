module OrdinaryDifferentialEquation

include("ordinary_differential_equation/test_methods.jl")
export test_method, test_pairs_method

include("ordinary_differential_equation/euler.jl")
export euler

include("ordinary_differential_equation/explicit_trap.jl")
export explicit_trap

include("ordinary_differential_equation/vectorization.jl")
export vectorization

include("ordinary_differential_equation/the_pendulum.jl")
export the_pendulum, the_undamped_pendulum_model, the_damped_pendulum_model, the_forced_damped_pendulum_model
export the_double_pendulum_model

include("ordinary_differential_equation/rk_midpoint.jl")
export rk_midpoint

include("ordinary_differential_equation/rk4.jl")
export rk4

include("ordinary_differential_equation/hodgkin_huxley.jl")
export hodgkin_huxley

end