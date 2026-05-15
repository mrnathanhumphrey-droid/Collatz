# IGUSA_C — Bories-Veys non-degenerated surface singularities

## Phase 0 — verbatim hypothesis

**Theorem 0.12 (Bories-Veys, Monodromy Conjecture for Igusa's p-adic local zeta function of a non-degenerated surface singularity).** Let f(x,y,z) ∈ Z[x,y,z] be a nonzero polynomial in **three variables** satisfying f(0,0,0) = 0. Suppose that f is non-degenerated over C with respect to all the compact faces of its Newton polyhedron, and let p be a prime such that f is also non-degenerated over F_p with respect to the same faces. If s_0 is a pole of Z_f^0, then exp(2πi Re(s_0)) is an eigenvalue of the local monodromy of f at some point of f^{-1}(0) ∩ U.

## Hypothesis types

- (i) Polynomial: **three variables**, f(0)=0, non-degenerate over C AND F_p w.r.t. compact faces.
- (ii) Conclusion: monodromy eigenvalues match exp(2πi Re(poles)).

## Phase 1 — substrate check

R78 substrate g(u) is **univariate** (n=1). Theorem requires **n=3**. **h_DIMENSION FAILED categorically.**

Even if we artificially lift to n=3 (e.g., by considering Z_3-valued g of three variables), the substrate has nothing to do with surface singularities.

## Disposition: NO_FIT (categorical, dimension mismatch)
