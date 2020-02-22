module OrdinaryDifferentialEquation

include("ordinary_differential_equation/euler.jl")
export euler, test_euler

include("ordinary_differential_equation/explicit_trap.jl")
export explicit_trap, test_explicit_trap

include("ordinary_differential_equation/vectorization.jl")
export vectorization

include("ordinary_differential_equation/the_pendulum.jl")
export the_pendulum, the_undamped_pendulum_model, the_damped_pendulum_model, the_forced_damped_pendulum_model
export the_double_pendulum_model

end