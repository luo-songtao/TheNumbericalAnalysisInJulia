include("naive_gauss_elimination.jl")
include("lu_factorization.jl")
include("lu_factorization_with_permutation.jl")
include("cholesky_decomposition.jl")
include("conjugate_gradient.jl")
include("conjugate_gradient_with_pre_condition.jl")
include("gauss_seidel.jl")
include("jacobi.jl")
include("sor.jl")
include("ssor.jl")

n = 10000
A = zeros(n,n)
A[1,1] = 3
for i = 2:n
    A[i,i] = 3
    A[i, i-1] = -1
    A[i-1, i] = -1
    
    if i != n/2 && i !=n/2+1
        A[i, n+1-i] = 0.5
    end
end

b = zeros(n)
for i = 2:n-1
    if i == n/2 || i ==n/2+1
        b[i] = 1.0
    else
        b[i] = 1.5
    end
end
b[1] = 2.5
b[n] = 2.5

# println(A)
# println(b)

x = @timev jacobi(A, b, zeros(n), n, n)
println(x[1])

# x = @timev sor(A, b, zeros(n), n, n, 1.1)
# println(x)

# x = @timev conjugate_gradient(A, b, zeros(n), n)
# println(x)

x = @timev ssor(A, b, zeros(n), n, 1)
println(x[1])