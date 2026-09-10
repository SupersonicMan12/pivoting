using LinearAlgebra

function complete_pivot_growth(A::AbstractMatrix)
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
    pivot_cols = Int[]
    
    U = Matrix{Float64}(A)
    # improved: just iterate through all 1-n, when k = n
    # it is effectively what we do specially at the end
    for k in 1:n
        maxr, maxc = k, k
        # tie-breaking in column-major, row-minor order
        for c in k:n, r in k:n
            if (maxr == -1 || abs(U[r, c]) > abs(U[maxr, maxc]))
                maxr, maxc = r, c
            end
        end
        stage_maximum = abs(U[maxr, maxc])
        if stage_maximum == 0.0
            throw(SingularException(k))
        end

        push!(stage_maxima, stage_maximum)
        push!(pivot_rows, maxr)
        push!(pivot_cols, maxc)

        if maxr != k
            U[[k, maxr], :] = U[[maxr, k], :]
        end
        if maxc != k
            U[:, [k, maxc]] = U[:, [maxc, k]]
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
        pivot_cols,
        max_multipliers,
        U
    )
end