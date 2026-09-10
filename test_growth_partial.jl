using LinearAlgebra
using Random
using Test
include("growth_partial.jl")
include("schur.jl")

@testset "Hand-written testcases" begin
    A = reshape([-3.0], 1, 1)
    result = partial_pivot_growth(A)
    @test result.initial_max == 3.0
    @test result.stage_maxima == [3.0]
    @test result.stage_growth == [1.0]
    @test result.growth == 1.0
    @test result.pivots == [-3.0]
    @test result.max_multipliers == [0.0]

    A = [
        1.0 1.0
        1.0 -1.0
    ]
    result = partial_pivot_growth(A)
    @test result.initial_max == 1.0
    @test result.stage_maxima == [1.0, 2.0]
    @test result.stage_growth == [1.0, 2.0]
    @test result.growth == 2.0
    @test result.pivots == [1.0, -2.0]
    @test result.max_multipliers == [1.0, 0.0]

    A = Matrix(Diagonal([4.0, 3.0, 2.0, 1.0]))
    result = partial_pivot_growth(A)
    @test result.initial_max == 4.0
    @test result.stage_maxima == [4.0, 3.0, 2.0, 1.0]
    @test result.stage_growth == [1.0, 0.75, 0.5, 0.25]
    @test result.growth == 1.0

    A = [
        1.0 2.0
        3.0 4.0
    ]

    # [
    # 3.0 4.0
    # 0 2.0-4.0/3.0
    # ]

    result = partial_pivot_growth(A)
    @test result.pivot_rows[1] == 2
    @test result.stage_maxima ≈ [4.0, 2.0/3.0]
    @test result.pivots ≈ [3.0, 2.0/3.0]
    @test result.max_multipliers ≈ [1.0/3.0, 0.0]
    @test result.growth ≈ 1.0
end

# same stuff as test_growth_complete.jl
@testset "Invalid testcases" begin
    @test_throws ArgumentError partial_pivot_growth(zeros(0, 0))
    @test_throws ArgumentError partial_pivot_growth(zeros(2, 3))
    @test_throws SingularException partial_pivot_growth(zeros(2, 2))
    @test_throws SingularException partial_pivot_growth([
        1.0 2.0
        2.0 4.0
    ])
end

# same stuff as test_growth_complete.jl
@testset "Random testcases: basic properties" begin
    rng = MersenneTwister(2009)
    for n in 1:25, trial in 1:20
        A = randn(rng, n, n)
        original = copy(A)
        result = partial_pivot_growth(A)
        @test A == original
        @test length(result.stage_maxima) == n
        @test length(result.stage_growth) == n
        @test length(result.pivots) == n
        @test length(result.pivot_rows) == n
        @test length(result.max_multipliers) == n
        @test result.initial_max == maximum(abs.(A))
        @test result.stage_maxima[1] == result.initial_max
        @test result.stage_growth ≈ result.stage_maxima ./ result.initial_max
        @test result.growth ≈ maximum(result.stage_growth)
        @test all(result.stage_maxima .> 0)
        @test all(result.stage_growth .> 0)
        @test result.growth >= 1
        @test all(abs.(result.pivots) .<= result.stage_maxima)
        @test all(x -> 0 <= x <= 1 + 16eps(Float64), 
            result.max_multipliers)
        @test all(k <= result.pivot_rows[k] <= n for k in 1:n)
        @test istriu(result.U)
    end
end

# same stuff as test_growth_complete.jl
@testset "Random testcases: scale invariant" begin
    rng = MersenneTwister(2010)
    for n in 1:15, trial in 1:10
        A = randn(rng, n, n)
        expected = partial_pivot_growth(A).growth
        @test partial_pivot_growth(0.5A).growth ≈ expected
        @test partial_pivot_growth(6.55A).growth ≈ expected
        @test partial_pivot_growth(-7.0A).growth ≈ expected
    end
end

@testset "Random testcases: Schur Cross-check" begin
    rng = MersenneTwister(2011)
    for n in 1:50, trial in 1:20
        A = randn(rng, n, n)
        expected = schur_partial_growth(A)
        result = partial_pivot_growth(A)
        @test result.growth ≈ expected.growth
        @test result.stage_maxima ≈ expected.stage_maxima
        @test result.stage_growth ≈ expected.stage_growth
    end
end
