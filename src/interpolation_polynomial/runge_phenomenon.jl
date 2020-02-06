"""
# 龙格例子

等距点插值会带来收敛困难当插值点数量增加。这一困难被称为龙格现象

"""
f(x) = 1/(1+12*x^2)

include("horner_rule.jl")
include("newton_difference_quotient.jl")
using Plots
pyplot()

nd = newton_difference_quotient
x_1 = -1:1/7:1
c_1 = nd(x_1, f.(x_1))
x_2 = -1:1/14:1
c_2 = nd(x_2, f.(x_2))

x = -1.2:0.1:1.2

s1 = scatter(x_1, f.(x_1), label="base point")
s2 = scatter(x_2, f.(x_2), label="base point")

p = plot(
    plot!(s1, x, [f.(x) horner_rule(length(x_1)-1, c_1, x, x_1)], label=["f" "f1"], ylims = (-1,2), yticks=-1:1:2, title="15 base point"),
    plot!(s2, x, [f.(x) horner_rule(length(x_2)-1, c_2, x, x_2)], label=["f" "f2"], ylims = (-1,2), yticks=-1:1:2, title="30 base point")
)
savefig(p, "a.png")