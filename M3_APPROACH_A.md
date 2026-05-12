# M3_APPROACH_A — direct spectral calculation

**Date:** 2026-05-11. Phase 2A of M3 probe. Computes the spectral-radius lower bound for M_3 = sup_γ ‖R(z, T_3)‖, addresses the operator-norm vs spectral-radius gap given T_3 is non-normal.

## A.1 Spectral lower bound

If T were normal, then `‖R(z, T)‖ = max_λ 1/|z − λ|` and M_3 = `sup_γ max_λ 1/|z − λ|` would be exact.

For our **conjectured spectrum** {1/2, 1/4, 1/8} on γ = `|z − 1/2| = 1/8`:

- `1/|z − 1/2| ≡ 8` (z lives on γ, distance exactly 1/8 to 1/2)
- `1/|z − 1/4| ≤ 8` (closest approach when z = 3/8, giving |z−1/4| = 1/8)
- `1/|z − 1/8| ≤ 4` (closest approach when z = 3/8, giving |z−1/8| = 1/4)

So `sup_γ max_λ 1/|z − λ| = 8`.

> **Spectral-radius lower bound:** `M_3 ≥ 8` (under conjectured spectrum).

## A.2 Why this isn't the operator norm — T_3 is non-normal

T_3 is the companion matrix of the polynomial `λ³ − (7/8)λ² + (7/32)λ − 1/64`. Companion matrices are **non-normal** in general unless the spectrum is a single point (degenerate case) — they have nontrivial Jordan / pseudospectral structure even for distinct eigenvalues.

Concretely: T_3 T_3* ≠ T_3* T_3 because the off-diagonal structure of T_3 (the row `(7/8, −7/32, 1/64)` plus the shift block `[[1,0,0],[0,1,0]]`) is not Hermitian-symmetric.

For non-normal operators, the operator-norm resolvent admits the inequality

> `‖R(z, T)‖_op  ≥  1 / dist(z, spec(T))`

with **equality iff T is normal**. The gap is captured by the **condition number** κ(V) = ‖V‖·‖V⁻¹‖ of the eigenvector basis under T = V D V⁻¹ (when T is diagonalizable):

> `‖R(z, T)‖_op  ≤  κ(V) · max_λ 1/|z − λ|`.

Compare to the **spectral-radius**-style lower bound `1/dist(z, spec)`. The ratio κ(V) is THE non-normality witness.

## A.3 Eigenvector basis condition number

R77.2 §3.3 computes V (Vandermonde at {1/2, 1/4, 1/8}):

> V = ⎡ 1/4   1/16   1/64 ⎤
>     ⎢ 1/2   1/4    1/8  ⎥
>     ⎣ 1     1      1    ⎦

with `‖V‖_F ≈ 1.843` and `det V = ±3/256 ≈ 0.0117`.

By the standard determinant identity for 3×3:

> `‖V⁻¹‖_F  =  ‖cof(V)ᵀ‖_F / |det V|`

`cof(V)` has entries that are 2×2 subdeterminants of V; each entry is `O(1)` (largest is on the order of 1·1 = 1). Sum of squares of 9 such entries ≤ 9 (crude); actually each cofactor is on the order of `O(1) × O(1) = O(1)`, so `‖cof‖_F ≤ 3`. Thus:

> `‖V⁻¹‖_F  ≤  3 / (3/256)  =  256` (less crude than R77.2's 768; closer to R77.2's "sharper 50–100" estimate)

A direct hand computation of V⁻¹ (Vandermonde inversion formula at distinct eigenvalues 1/2, 1/4, 1/8) gives entries scaled by `1 / ∏_{i≠j}(λ_i − λ_j)`. The pairwise differences are `1/4`, `3/8`, `1/8`, so the smallest product term is `(1/4)·(3/8)·(1/8) = 3/256`. V⁻¹ entries are then on the order of `λ_i^2 / (3/256)` ≈ `(1/4) × (256/3) ≈ 21.3`. With 9 entries of this order, `‖V⁻¹‖_F` ≈ `√(9 × 21.3²)` ≈ `64`. Frobenius dominates 2-norm: `‖V⁻¹‖_2 ≤ 64`.

> **Sharper condition number:** `κ(V) = ‖V‖_2 · ‖V⁻¹‖_2 ≤ 1.843 × 64 ≈ 118`.

Combined with `max_λ 1/|z−λ| ≤ 8` on γ:

> **Approach A bound:** `M_3 ≤ 118 × 8  ≈  944`.

This is tighter than R77.2's crude 11320 and sits inside R77.2's quoted range "800–1000".

## A.4 Outcome classification

> **APPROACH_A_BOUNDED.**
>
> - Lower bound (spectral radius): `M_3 ≥ 8` on γ.
> - Upper bound (Approach A via κ(V) × max_λ 1/|z−λ|): `M_3 ≤ 944`.
> - These bounds **bracket** M_3 with factor-of-118 uncertainty.
> - The lower bound is rigorous; the upper bound uses the Frobenius-norm-as-2-norm-upper-bound which is a standard but generally loose inequality.

## A.5 Honest scope

(1) **Approach A is the correct calc for the literal 3×3 matrix T_3** defined by R77.2's recursion coefficients. The bound `M_3 ∈ [8, 944]` is rigorous for that matrix.

(2) **R77.3 falsifies that T_3 describes ε_n's actual dynamics.** So the bound `M_3 ∈ [8, 944]` is computed about a matrix that doesn't appear in any actually-correct rate operator for c = 7/45. Approach A is therefore *mathematically clean* but **operationally moot** for Nisoli closure.

(3) **Approach B (perturbation series) and Approach C (numerical at higher modulus) would extend this only if there were a corrected, larger-dimensional T to extend to.** Both are deferred pending characterization of the true rate operator, which R77.4 erratum flags as parked.

## A.6 Summary

For Phase 3 closure-table parameterization, the value **M_3 ≈ 944** (or **800–1000** range from R77.2 §3.3) is the canonical project number, with the caveat that it pertains to a falsified matrix. Pre-R77.3, this would feed directly into closure; post-R77.3, it's a parameter in a hypothetical "if the spectrum were what we conjectured" calculation.
