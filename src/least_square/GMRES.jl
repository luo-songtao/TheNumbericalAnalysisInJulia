"""
# 广义最小余项方法(GMRES)
    gmres(A, b, x_0)
共轭梯度法可以看作一种迭代方法用于求解A为对称方阵的方程组Ax=b的解，但如果A不对称，则不能使用共轭梯度法。

但GMRES方法可用于求解非对称的矩阵A，而且它是求解大规模、稀疏、非对称线性方程组Ax=b的好方法。

在GMRES方法中涉及正交性来进行求解：
- 在每个迭代过程中使用最小二乘公式，最小化后向误差
- 在每步中搜索空间的基被重新正交化以消除病态问题中的不精确

GMRES属于Krylov方法，它依赖精确的Krylov空间计算，该空间是向量``\\{ r, Ar, ..., A^kr \\}``所张成的空间，，其中``r=b-Ax_0``。

由于向量``A^kr``对于大的k倾向一个共有方向，Krylov空间的基必须认真计算，找出Krylov空间精确的基需要正交化方法，例如格拉姆-施密特正交或这豪斯霍尔德反射子方法。
(本实现是基于格拉姆-施密特正交方法)

迭代过程相关公式说明：
- 矩阵``A``为``n\\times n``、 ``H``为``k+1\\times k``、``Q_k``为``n \\times k``，大多数情况下，k比n小很多
- ``x_k = x_add + x_0``, ``x_0``为初识估计，``x_k``为第k步的预估值，``x_{add}``为从``Q_k``中搜索到的用于改进原始估计的量
- ``x_add = Q_kc``
- ``AQ_k = Q_{k+1}H_k``

``\\begin{aligned} \\qquad 为了最小化Ax&=b的余项r \\\\ \\Vert b-A(x_0+x_{add}) \\Vert &= \\Vert r-Ax_{add} \\Vert \\\\ &= \\Vert AQ_kc - r \\Vert &= \\Vert Q_{k+1}H_kc - r \\Vert &= \\Vert H_kc-Q^T_{k+1}r \\Vert \\end{aligned}``

因此算法中后部分会利用最小二乘求解``c_k``，然后求出``x_add``

## 重启GMRES

GMRES方法的变种-重启GMRES，其中一种算法实现思想：如果第k步迭代后没能足够趋近解，
而且``n\\times k 的Q_k``变得大得难以处理，可以考虑扔掉``Q_k``，重新开始GMRES方法，并使用当前的最优估计``x_k``作为新的``x_0``

# Example
```jldoctest
julia> gmres([1 1 0; -1 1 2; 0 0 1], [1;0;0], [0.0;0.0;0.0])
3-element Array{Float64,1}:
 0.5
 0.5
 0.0
julia> gmres([1 1 0; 0 1 0; 1 1 1], [1;0;0], [0.0;0.0;0.0])
3-element Array{Float64,1}:
  1.0
  0.0
 -1.0
julia> gmres([0 0 1; 1 0 0; 0 1 0], [1;0;0], [0.0;0.0;0.0])
3-element Array{Float64,1}:
 0.0
 0.0
 1.0
```
"""
function gmres(A, b, x_0)
    m,n = size(A)
    r = b - A*x_0
    q = r/ sqrt(sum(r.^2))

    H = zeros(m+1, m)
    Q = zeros(m, m+1)
    Q[:, 1] = q
    x_k = [x_0]
    for k = 1:m
        y = A*Q[:, k]
        for j = 1:k
            H[j,k] = Q[:, j]'*y
            y = y - H[j,k] * Q[:, j]
        end
        H[k+1, k] = sqrt(sum(y.^2))
        if H[k+1, k] != 0
            Q[:, k+1] = y / H[k+1, k]
        end
        c_0 = zeros(k); br= zeros(k+1); h = zeros(k+1)
        br[1] = sqrt(sum(r.^2)); h[1] = 1.0
        c_k = least_square_by_complete_qr(hcat(H[1:k+1, 1:k], h), br, c_0)
        x_k[1] = Q[:, 1:k]*c_k + x_0
        if H[k+1, k] == 0
            break
        end
    end
    return x_k[1]
end
             
push!(LOAD_PATH, "../")
using LeastSquare

"""
# 预条件GMRES
    pre_con_gmres(A, b, x_0, M)

同共轭梯度法类似，使用：

```math
M^{-1}Ax = M^{-1}b
```

``M``被称为预条件子

# Example
```jldoctest
julia> pre_con_gmres([1 1 0; -1 1 2; 0 0 1], [1;0;0], [0.0;0.0;0.0], [1 0 0;0 1 0; 0 0 1])
3-element Array{Float64,1}:
 0.5
 0.5
 0.0
```
"""
function pre_con_gmres(A, b, x_0, M)
    m,n = size(A)
    M_inv = inv(M)
    r = M_inv*(b - A*x_0)
    q = r/ sqrt(sum(r.^2))

    H = zeros(m+1, m)
    Q = zeros(m, m+1)
    Q[:, 1] = q
    x_k = [x_0]
    for k = 1:m
        y = M_inv*A*Q[:, k]
        for j = 1:k
            H[j,k] = Q[:, j]'*y
            y = y - H[j,k] * Q[:, j]
        end
        H[k+1, k] = sqrt(sum(y.^2))
        if H[k+1, k] != 0
            Q[:, k+1] = y / H[k+1, k]
        end
        c_0 = zeros(k); br= zeros(k+1); h = zeros(k+1)
        br[1] = sqrt(sum(r.^2)); h[1] = 1.0
        c_k = least_square_by_complete_qr(hcat(H[1:k+1, 1:k], h), br, c_0)
        x_k[1] = Q[:, 1:k]*c_k + x_0
        if H[k+1, k] == 0
            break
        end
    end
    return x_k[1]
end
