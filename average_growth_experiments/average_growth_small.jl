# overall route:
# julia: generate calculate save to csv
# csv as the bridge
# python does data processing + visualizing

# run:
# julia average_growth_small.jl
# python3 analysis.py

using Random
include(joinpath(@__DIR__, "..", "growth_partial.jl"))
include(joinpath(@__DIR__, "..", "growth_complete.jl"))
sizes = [16, 32, 64]
trials = 20
base_seed = 2009
output_path = joinpath(@__DIR__, "pilot_growth.csv")

open(output_path, "w") do io
    println(io, "n,trial,matrix_seed,partial_growth,complete_growth")
    for n in sizes, trial in 1:trials
        matrix_seed = base_seed + 1000000n + trial
        rng = MersenneTwister(matrix_seed)
        A = randn(rng, n, n)
        partial = partial_pivot_growth(A).growth
        complete = complete_pivot_growth(A).growth
        println(io,"$n,$trial,$matrix_seed,$partial,$complete")
    end
end
