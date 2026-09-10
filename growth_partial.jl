using LinearAlgebra

function partial_pivot_growth(A::AbstractMatrix)
    n, m = size(A)
    if n != m || n == 0 || m == 0
        throw(ArgumentError("A must be square"))
    end

    initial_max = maximum(abs.(A))
    if initial_max == 0.0
        throw(SingularException(0))
    end

    # necessary:
    stage_maxima = Float64[]
    # additional stats for reproducibility 
    pivots = Float64[]
    max_multipliers = Float64[]
    pivot_rows = Int[]
    
    U = Matrix{Float64}(A)
    for k in 1:n
        maxr = k
        for r in k:n
            if (abs(U[r, k]) > abs(U[maxr, k]))
                maxr = r
            end
        end
        pivot_magnitude = abs(U[maxr, k])
        if pivot_magnitude == 0.0
            throw(SingularException(k))
        end

        stage_maximum = 0.0
        for c in k:n, r in k:n
            stage_maximum = max(stage_maximum, abs(U[r, c]))
        end
        
        push!(stage_maxima, stage_maximum)
        push!(pivot_rows, maxr)

        if maxr != k
            U[[k, maxr], :] = U[[maxr, k], :]
        end

        push!(pivots, U[k, k])
        largest_multiplier = 0.0

        for i in k+1:n
            multiplier = U[i,k] / U[k,k]
            largest_multiplier = max(largest_multiplier, abs(multiplier))
            for j in k:n
                U[i, j] -= multiplier * U[k, j]
            end
            U[i, k] = 0.0
        end
        push!(max_multipliers, largest_multiplier)
    end

    stage_growth = stage_maxima ./ initial_max;
    growth = maximum(stage_growth)

    return ( ;
        growth,
        initial_max,
        stage_maxima,
        stage_growth,
        pivots,
        pivot_rows,
        max_multipliers,
        U
    )
end