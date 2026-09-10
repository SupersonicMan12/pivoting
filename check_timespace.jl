include("growth_complete.jl")
using Random

rng = MersenneTwister(2009)

for n in (64, 128, 256, 512)
    A = randn(rng, n, n)
    stats = @timed complete_pivot_growth(A)

    println(
        "n=$n ",
        "time=$(round(stats.time, digits=4))s ",
        "memory=$(round(stats.bytes / 2.0^20, digits=2)) MiB ",
        "growth=$(round(stats.value.growth, digits=4))",
    )
end