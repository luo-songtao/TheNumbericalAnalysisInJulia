module LeastSquare

include("least_square/normal_equation.jl")
export normal_equation, fit_modeling_by_normal_equation

include("least_square/gram_schmidt_orthogon.jl")
export classic_gram_schmidt_orthogon, gram_schmidt_orthogon

include("least_square/least_square_and_qr.jl")
export incomplete_qr,  complete_qr, least_square_by_complete_qr

include("least_square/householder_reflector.jl")
export householder_reflector, qr_by_householder_reflector

end