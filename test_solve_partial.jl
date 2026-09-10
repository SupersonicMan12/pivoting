using LinearAlgebra
using Random
using Test
include("solve_partial.jl")

function check_solve(A::AbstractMatrix, rng::AbstractRNG)
    n = size(A, 1)
    fact = initialize(A)
    @test fact.P * A ≈ fact.L * fact.U
    b = randn(rng, n)
    x = solve(fact, b)

    # Residual check: ||Ax - b|| has relative error & abs error within 1e-12
    residual = norm(A*x-b)
    scale = max(norm(A)*norm(x), norm(b), 1)
    @test residual <= 1e-12 * scale
end

@testset "Ax=b Solve" begin
    @testset "Randomized testcases" begin
        rng = MersenneTwister(2009)
        for n in 1:100, trial in 1:20
            @testset "n = $n, trial = $trial" begin
                check_solve(randn(rng, n, n), rng)
            end
        end
    end
end
