# overall route:
# julia: generate calculate save to csv
# csv as the bridge
# python does data processing + visualizing

# run:
# julia average_growth.jl
# python3 analysis.py

# runtime roughly 3-4 minutes!!!!
# total expected: 7240 test cases

using Random
include(joinpath(@__DIR__, "..", "growth_partial.jl"))
include(joinpath(@__DIR__, "..", "growth_complete.jl"))
sizes = [16, 32, 64, 128, 256, 512, 1024]
trial_by_size = Dict(
    16 => 2000,
    32 => 2000,
    64 => 1500,
    128 => 1000,
    256 => 500,
    512 => 200,
    1024 => 40
)
base_seed = 2009
output_path = joinpath(@__DIR__, "gaussian_growth.csv")

open(output_path, "w") do io
    println(io, "n,trial,matrix_seed,partial_growth,complete_growth")
    for n in sizes
        trials = trial_by_size[n]
        for trial in 1:trials
            matrix_seed = base_seed + 1000000n + trial
            rng = MersenneTwister(matrix_seed)
            A = randn(rng, n, n)
            partial = partial_pivot_growth(A).growth
            complete = complete_pivot_growth(A).growth
            println(
                io,
                "$n,$trial,$matrix_seed,$partial,$complete"
            )
        end
        flush(io)
        println("Finished n = $n with $trials trials")
    end
end
