# use this to check complete pivoting implementation

using LinearAlgebra

function schur_complete_growth(A::AbstractMatrix)
    S = Matrix{Float64}(A)
    initial_max = maximum(abs.(S))
    stage_maxima = Float64[]
    while size(S, 1) > 0
        pivot_index = argmax(abs.(S))
        pivot_row, pivot_col = Tuple(pivot_index)
        stage_maximum = abs(S[pivot_row, pivot_col])
        if stage_maximum == 0.0
            throw(SingularException(length(stage_maxima) + 1))
        end
        push!(stage_maxima, stage_maximum)
        if pivot_row != 1
            S[[1, pivot_row], :] = S[[pivot_row, 1], :]
        end
        if pivot_col != 1
            S[:, [1, pivot_col]] = S[:, [pivot_col, 1]]
        end
        if size(S, 1) == 1
            break
        end
        pivot = S[1, 1]
        column = S[2:end, 1]
        row = S[1, 2:end]
        B = S[2:end, 2:end]
        # suppose [pivot row; column B]
        # S new = B - column * row transpose / pivot 
        S = B - column * transpose(row) / pivot
    end
    stage_growth = stage_maxima ./ initial_max
    growth = maximum(stage_growth)
    return (;
        growth, 
        stage_maxima,
        stage_growth
    )
end

function schur_partial_growth(A::AbstractMatrix)
    S = Matrix{Float64}(A)
    initial_max = maximum(abs.(S))
    stage_maxima = Float64[]
    while size(S, 1) > 0
        stage_maximum = maximum(abs.(S))
        push!(stage_maxima, stage_maximum)
        pivot_row = argmax(abs.(S[:, 1]))
        pivot_magnitude = abs(S[pivot_row, 1])
        if pivot_magnitude == 0.0
            throw(SingularException(length(stage_maxima)))
        end
        if pivot_row != 1
            S[[1, pivot_row], :] = S[[pivot_row, 1], :]
        end
        if size(S, 1) == 1
            break
        end
        pivot = S[1, 1]
        column = S[2:end, 1]
        row = S[1, 2:end]
        B = S[2:end, 2:end]
        S = B - column * transpose(row) / pivot
    end
    stage_growth = stage_maxima ./ initial_max
    growth = maximum(stage_growth)
    return (;
        growth, 
        stage_maxima,
        stage_growth
    )
end