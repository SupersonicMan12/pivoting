# Gaussian Elimination and Growth Factors

Implementations of Gaussian elimination with no pivoting, partial pivoting and complete pivoting in Julia.

The complete-pivoting growth tracker records the maximum absolute entry at every elimination stage and computes

\[
\rho(A) =
\frac{\max_k \|A^{(k)}\|_{\max}}
     {\|A\|_{\max}}.
\]

It is independently checked against various invariants and a separate Schur-complement implementation.

Tests:

```bash
julia --startup-file=no test_complete.jl
julia --startup-file=no test_partial.jl
julia --startup-file=no test_solve_partial.jl
julia --startup-file=no test_growth_complete.jl
julia --startup-file=no test_growth_partial.jl
```

## Average-growth experiment

Generate the paired Gaussian-matrix data:

```bash
julia --startup-file=no experiments/average_growth.jl
```

Validate and generate the median-growth figure:

```bash
python3 experiments/graph.py
```

This produced:

```markdown
![Median growth under partial and complete pivoting](figures/median_growth.png)
```

## Development Note

I implemented the elimination and growth-tracking routines while learning numerical linear algebra and Julia. I used an AI coding assistant, GPT 5.6 Sol Medium, for syntax guidance, code review, test and experiment planning and a plotting scaffold. All numerical results were generated locally. All implementations, except visualization, are completely handwritten.