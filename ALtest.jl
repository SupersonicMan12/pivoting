using LinearAlgebra
include("elimination.jl")

A = [
     2.0   1.0  -1.0
    -3.0  -1.0   2.0
    -2.0   1.0   2.0
]

L, U = lu_no_pivot(A)

display(L)
display(U)

println("Reconstruction error: ", norm(A - L * U))
println("Correct: ", A ≈ L * U)