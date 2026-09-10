using LinearAlgebra

# assume that no swaps are necessary
function verify_prescribed_complete_growth(A::AbstractMatrix)
    n, m = size(A)
    if n != m || n == 0 || m == 0
        throw(ArgumentError("A must be nonempty and square"))
    end

    S = Matrix(A)
    T = eltype(S)
    if !(T <: Rational)
        throw(ArgumentError("A must be at least rational"))
    end

    initial_max = maximum(abs.(S))
    stage_maxima = T[]
    pivots = T[]

    while size(S, 1) > 0
        pivot = S[1, 1]
        active_maximum = maximum(abs.(S))
        stage = length(pivots) + 1
        if iszero(pivot)
            throw(SingularException(stage))
        end
        if abs(pivot) != active_maximum
            throw(ArgumentError(
                "Prescribed pivot at stage $stage is not maximal as assumed"  
            ))
        end
        push!(pivots, pivot)
        push!(stage_maxima, active_maximum)
        if size(S, 1) == 1
            break
        end
        column = S[2:end, 1]
        row = S[1, 2:end]
        B = S[2:end, 2:end]
        S = B - column*transpose(row)/pivot
    end
    stage_growth = stage_maxima ./ initial_max
    growth = maximum(stage_growth)
    return (;
        growth,
        initial_max,
        stage_maxima,
        stage_growth,
        pivots
    )
end
    