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
julia --startup-file=no test_growth_complete.jl
julia --startup-file=no test_complete.jl
julia --startup-file=no test_partial.jl
julia --startup-file=no test_partial_solve.jl
```
