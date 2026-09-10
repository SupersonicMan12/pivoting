using LinearAlgebra
using Random
include("elimination.jl")

struct LUPFactorization{T<:Real}
    L::Matrix{T}
    U::Matrix{T}
    P::Matrix{T}
end

function initialize(A::AbstractMatrix)
    L, U, P = lup_partial_pivot(A)
    return LUPFactorization(L, U, P)
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

# We have PA = LU and Ax = b
# PAx = Pb
# LUx = Pb (1)
# let y = Ux, then y is the vector that has Ly = (Pb) (2)
# x is the vector that has Ux = y (3)
function solve(fact::LUPFactorization, b::AbstractVector)
    n = size(fact.L, 1)
    if length(b) != n
        throw(DimensionMismatch("Vector length $(length(b)) does not match matrix size $n"))
    end
    c = fact.P * b
    y = forward_substitution(fact.L, c)
    x = back_substitution(fact.U, y)
    return x
end