using JLD

include("exact_growth.jl")

path = joinpath(@__DIR__, "n05_4p1325.jld")
data = load(path)
B = data["B"]

result = verify_prescribed_complete_growth(B)

println("Exact rational growth:")
println(result.growth)
println()
println("Decimal growth:")
println(BigFloat(result.growth))
println()
println("Stage maxima")
for (stage, maximum_value) in enumerate(result.stage_maxima)
    println("stage $stage: ",
    BigFloat(maximum_value))
end
println()
println("Pivots:")
for (stage, pivot) in enumerate(result.pivots)
    println("stage $stage: ",
    BigFloat(pivot))
end
