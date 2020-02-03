"""
# 对称连续过松弛(SSOR)预条件共轭梯度法
    ssor(A, b, x_0, n, ω=1)

``M = (D + \\omega L)D^{-1} (D + \\omega U)``

同SOR方法中，``A = L+D+U``，``\\omega``是0和2之间的常数。

当``\\omega =1 ``时，它也被称为**高斯-塞德尔预条件子**

``\\begin{aligned} M &= (D + \\omega L)D^{-1} (D + \\omega U) \\\\ M &= (I + \\omega LD^{-1}) (D + \\omega U) \\\\ z &= M^{-1}v \\\\ Mz &= v \\\\ & (I + \\omega LD^{-1})c = v \\\\ & (D + \\omega U)z = c \\end{aligned}``

对于稀疏矩阵，两次回代所花的时间和非零元素的个数成正比
"""
function ssor(A, b, x_0, n, ω=1)
    x_k = x_0
    r_k = b - A*x_k    # 这里是r_0
    z_k = update_z_vertor(A, r_k, n, ω)    # 这里是z_0
    d_k = z_k    # 这里是d_0
    for k = 0:n-1
        if iszero(r_k); break; end
        α_k = (transpose(r_k) * z_k)/(transpose(d_k)*A*d_k)
        x_k = x_k + α_k * d_k
        r_k_plus_1 = r_k - α_k* A * d_k
        z_k_plus_1 = update_z_vertor(A, r_k_plus_1, n, ω)
        β_k = (transpose(r_k_plus_1)*z_k_plus_1)/(transpose(r_k)*z_k)
        d_k = z_k_plus_1 + β_k * d_k
        r_k = r_k_plus_1
        z_k = z_k_plus_1
    end
    return x_k
end

function update_z_vertor(A, r, n, ω)
    v = deepcopy(r)
    c = zeros(n)
    # 第一次回代解出c
    for i = 1:n
        for j = 1:i-1
            v[i] = v[i] - (ω * A[i,j] * 1/A[j,j])*c[j]
        end
        c[i] = v[i]
    end
    z = zeros(n)
    # 第二次回代解出z
    for i = n:-1:1
        for j = i+1:n
            c[i] = c[i] - (ω * A[i,j])*z[j]
        end
        z[i] = c[i]/A[i,i]
    end
    return z
end

