# Verifying the Prescribed Matrix that Achieves 4.1325170311...

This directory independently evaluates the rational $5\times5$ witness published in the companion repository to Edelman and Urschel's work on complete-pivoting growth (https://github.com/alanedelman/CompletePivotingGrowth/tree/main)

The source file is [`Matrices/n05_4p1325.jld`](https://github.com/alanedelman/CompletePivotingGrowth/blob/main/Matrices/n05_4p1325.jld)

The JLD file contains two matrices:
`A`, a floating-point numerical optimizer output
`B`, an exact witness derived from `A` and repaired so that its prescribed diagonal pivots are maximal at every stage.

This verification uses `B`.

## Method

The matrix B already has the intended row and column permutation, hence no swaps are needed. At each elimination step, we

1. take the top left entry as the pivot
2. check if the pivot is maximal in the active matrix
3. form the next Schur complement
4. record the exact pivot and stage maximum

All computations are done in `Rational{BigInt}`

## Result

All five prescribed pivots of `B` are legal, and the exact growth factor is a rational number whose decimal value is

```output
Decimal growth:
4.132517031148570089067536169627026009232245454614513962084506271889200206000998
```

This verifies the lower bound of $$g_5 \geq 4.132517031148570089\ldots$$

The best currently established upper bound is approximately $4.84$ (Chen, Edelman and Urschel 2026). 

## Run

```bash
julia --startup-file=no --project=g5 g5/verify_n5.jl
```