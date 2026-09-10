# no pivoting / A = LU
# we assume what we get is square, non-singular and decomposable

using LinearAlgebra

function lu_no_pivot(A::AbstractMatrix)
    n, m = size(A)
    if n != m
        throw(ArgumentError("A must be square"))
    end
    U = Matrix{Float64}(A)
    L = Matrix{Float64}(I, n, n)
    for k in 1:n-1
        if U[k,k] == 0.0
            # Not necessarily singular!
            error("Unexpected Zero Pivot")
        end
        for i in k+1:n
            multiplier = U[i,k]/U[k,k]
            L[i,k] = multiplier
            for j in k:n
                U[i, j] -= multiplier * U[k, j]
            end
            # force 0.0
            U[i, k] = 0.0
        end
    end
    if U[n, n] == 0.0
        throw(SingularException(n))
    end
    return L, U
end

# (QP)(A) = (QLQ)(QU)
# not using permutation matrix just yet
# "search time" = how long total to find partial pivot 
# O(n^2)
function lup_partial_pivot(A::AbstractMatrix)
    n, m = size(A)
    if n != m
        throw(ArgumentError("A must be square"))
    end
    U = Matrix{Float64}(A)
    L = Matrix{Float64}(I, n, n)
    P = Matrix{Float64}(I, n, n)
    for k in 1:n-1
        maxr = k-1+argmax(abs.(U[k:n, k]))
        if U[maxr, k] == 0
            throw(SingularException(k))
        end
        # swap rows: U swap, P swap,
        # L experiences transpose on the subset of 
        # row maxr, k union column maxr, k
        if k != maxr
            U[[k, maxr], :] = U[[maxr, k], :]
            P[[k, maxr], :] = P[[maxr, k], :]
            
            # QLQ^T = QLQ yields the following
            # technically:
            # L[[k, maxr], :] = L[[maxr, k], :]
            # L[:, [k, maxr]] = L[:, [maxr, k]]

            # this is mathematically equivalent to:
            L[[k, maxr], 1:k-1] = L[[maxr, k], 1:k-1]
        end
        for i in k+1:n
            multiplier = U[i,k] / U[k,k]
            L[i, k] = multiplier
            for j in k:n
                U[i, j] -= multiplier * U[k, j]
            end
            U[i, k] = 0.0
        end
    end
    if U[n, n] == 0.0
        throw(SingularException(n))
    end
    return L, U, P
end

# row swap: (SP)(A)Q = (SLS)(SU)
# column swap: PA(QS) = L(US)
# search time O(n^3)
function lupq_complete_pivot(A::AbstractMatrix)
    n, m = size(A)
    if n != m
        throw(ArgumentError("A must be square"))
    end
    U = Matrix{Float64}(A)
    L = Matrix{Float64}(I, n, n)
    P = Matrix{Float64}(I, n, n)
    Q = Matrix{Float64}(I, n, n)
    for k in 1:n-1
        maxr, maxc = k, k
        for r in k:n, c in k:n
            if (maxr == -1 || abs(U[r, c]) > abs(U[maxr, maxc]))
                maxr, maxc = r, c
            end
        end
        if U[maxr, maxc] == 0.0
            throw(SingularException(k))
        end
        if maxr != k
            U[[k, maxr], :] = U[[maxr, k], :]
            P[[k, maxr], :] = P[[maxr, k], :]
            L[[k, maxr], 1:k-1] = L[[maxr, k], 1:k-1]
        end
        if maxc != k
            U[:, [k, maxc]] = U[:, [maxc, k]]
            Q[:, [k, maxc]] = Q[:, [maxc, k]]
        end
        for i in k+1:n
            multiplier = U[i,k] / U[k,k]
            L[i, k] = multiplier
            for j in k:n
                U[i, j] -= multiplier * U[k, j]
            end
            U[i, k] = 0.0
        end
    end
    if U[n, n] == 0.0
        throw(SingularException(n))
    end
    return L, U, P, Q
end
