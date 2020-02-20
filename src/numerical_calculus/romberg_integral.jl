"""
# Romberg Integral
    romberg_integral(f, a, b, n)
Romberg Integral是对复合梯形法则应用外推的结果。复合梯形法则是关于h的二阶法则，使用外推至少可以得到三阶的法则。

# Example
```jldoctest
julia> romberg_integral(x->log(exp(1), x), 1, 2, 4)
4×4 Array{Float64,2}:
 0.346574  0.0       0.0       0.0     
 0.376019  0.385835  0.0       0.0     
 0.3837    0.38626   0.386288  0.0     
 0.385644  0.386292  0.386294  0.386294
```
"""
function romberg_integral(f, a, b, n)
    h = b - a
    R = zeros(n,n)
    R[1,1] = h*(f(a)+f(b))/2
    for j = 2:n
        h = h/2
        i = 1:2^(j-2)
        R[j,1] = (1/2)*R[j-1,1] + h * sum(f.(h*(2i.-1).+a))
        for k = 2:j
            R[j,k] = (4^(k-1)*R[j,k-1] - R[j-1,k-1])/(4^(k-1)-1)
        end
    end
    return R
end

