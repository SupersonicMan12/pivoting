# Testing Script manually written,
# guided by GPT Astra 6.0 for the purpose of learning

using LinearAlgebra
using Test
using Random
include("elimination.jl")

function check_lup(A)
    original = copy(A)
    L, U, P = lup_partial_pivot(A)
    n = size(A, 1)
    @test A == original
    @test size(L) == size(U) == size(P) == size(A)
    @test all(isfinite, L) && all(isfinite, U) && all(isfinite, P)
    @test istril(L)
    @test istriu(U)
    @test diag(L) == ones(n)
    @test all(x -> x == 0 || x == 1, P)
    # dims = 1 over rows, dims = 2 over columns
    @test vec(sum(P; dims = 1)) == ones(n)
    @test vec(sum(P; dims = 2)) == ones(n)
    # tril(A, k) = k diagonals above the main diagonal, and downwards
    @test maximum(abs.(tril(L, -1))) <= 1+1e-12
    
    @test isapprox(P*A, L*U; rtol=1e-12, atol=0)
end

@testset "Partial-pivot LU" begin
    @testset "Targeted testcases" begin
        check_lup(reshape([3.0], 1, 1))
        check_lup(Matrix{Float64}(I, 2, 2))
        check_lup(Matrix{Float64}(I, 3, 3))
        check_lup([0 1;2 3])
        check_lup([0 1 0; 0 0 1; 1 0 0])
        check_lup([4 1 0; 2 1 1; 1 3 1])
        check_lup([1 2; -1 2])
        for scale in (1e-100, 1.0, 1e100)
            # automatic promotion of all the other elements
            # thanks to my 4.0
            check_lup(scale .* [4.0 1 0; 2 1 1; 1 3 1])
        end
    end

    @testset "Random matrices" begin
        rng = MersenneTwister(2009)
        for n in (1, 2, 3, 5, 10, 25, 50), trial in 1:20
            @testset "n = $n trial = $trial" begin
                check_lup(randn(rng, n, n))
            end
        end
    end

    @testset "Invalid Inputs" begin
        @test_throws ArgumentError lup_partial_pivot(zeros(2,3))
        for A in (
            zeros(1, 1),
            [0.0 1; 0 2],
            [1.0 2; 2 4],
            Matrix(Diagonal([1.0, 0, 1])),
        )
            @test_throws ErrorException lup_partial_pivot(A)
        end
    end
end