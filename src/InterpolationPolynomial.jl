module InterpolationPolynomial

include("interpolation_polynomial/horner_rule.jl")
export horner_rule

include("interpolation_polynomial/newton_difference_quotient.jl")
export newton_difference_quotient

include("interpolation_polynomial/chebyshev_interpolation.jl")
export chebyshev_interpolation

include("interpolation_polynomial/natural_cubic_spline.jl")
export natural_cubic_spline

include("interpolation_polynomial/curvature_adjustable_cubic_spline.jl")
export curvature_adjustable_cubic_spline

include("interpolation_polynomial/clamp_cubic_spline.jl")
export clamp_cubic_spline

include("interpolation_polynomial/parabola_endpoint_cubic_spline.jl")
export parabola_endpoint_cubic_spline

include("interpolation_polynomial/unknot_cubic_spline.jl")
export unknot_cubic_spline

include("interpolation_polynomial/bezier_curve.jl")
export bezier_curve, bezier_draw

end