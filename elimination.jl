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
            U[i, k] = 0.0
        end
    end
    if U[n, n] == 0.0
        error("Singular Matrix")
    end
    return L, U
end

# (QP)(A) = (QLQ)(QU)
# not using permutation matrix just yet
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
            error("Singular Matrix")
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
        error("Singular Matrix")
    end
    return L, U, P
end