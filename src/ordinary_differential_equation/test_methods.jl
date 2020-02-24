function test_method(method)
    y(t) = 3*exp(1)^(t^2/2)-t^2-2
    dy(t, y) = t*y + t^3
    a,b = [0,1]
    result = zeros(10,3)
    n = 5
    for i = 1:10
        h = 1/n
        err = abs(y(1) - method(dy, a, b, 1, n))
        result[i, :] = [n h err]
        n = 2n
    end
    return result
end
