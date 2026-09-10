# Analysis
# PAQ = LU, Ax = b
# PA = LUQ_T
# PAx = Pb
# LUQ_T x = Pb
# Compute Pb
# Forward substitution find UQ_T x
# Back substitution find Q_T x = Y
# then x = QY after multiplying

using LinearAlgebra
using Random
include("elimination.jl")

struct LUPQFactorization{T<:Real}
    L::Matrix{T}
    U::Matrix{T}
    P::Matrix{T}
    Q::Matrix{T}
end

function initialize(A::AbstractMatrix)
    L, U, P, Q = lupq_complete_pivot(A)
    return LUPQFactorization(L, U, P, Q)
end

# Lx = b
function forward_substitution(L::AbstractMatrix, b::AbstractVector)
    n = length(b)
    y = zeros(Float64, n)
    for i in 1:n
        y[i] = b[i] - dot(L[i, 1:i-1], y[1:i-1])
    end
    return y
end

# Ux = b
function back_substitution(U::AbstractMatrix, b::AbstractVector)
    n = length(b)
    y = zeros(Float64, n)
    for i in n:-1:1
        y[i] = (b[i] - dot(U[i, i+1:n], y[i+1:n])) / U[i, i]
    end
    return y
end 

function solve(fact::LUPQFactorization, b::AbstractVector)
    n = size(fact.L, 1)
    if length(b) != n
        throw(DimensionMismatch("Vector length $(length(b)) does not match matrix size $n"))
    end
    c = fact.P * b
    y = forward_substitution(fact.L, c)
    qtx = back_substitution(fact.U, y)
    x = fact.Q * qtx
    return x
end


