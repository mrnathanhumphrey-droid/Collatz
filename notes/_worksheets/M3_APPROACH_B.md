# M3_APPROACH_B — perturbation series

**Date:** 2026-05-11. Phase 2B of M3 probe. Tests whether T_3 admits a useful Neumann-series decomposition T_3 = T_0 + ε·T_1 with explicitly bounded resolvent.

## B.1 Setup

Approach B would write `T_3 = T_0 + ε · T_1` where:

- `T_0` is a chosen "base" operator with explicitly computable resolvent,
- `T_1` is the perturbation,
- ε is a small parameter (often 1, with smallness handled via `T_1` norm).

Then for `z ∈ γ`:

> `R(z, T_3) = R(z, T_0) · [I − ε · T_1 · R(z, T_0)]⁻¹  =  R(z, T_0) · Σ_{n≥0} (ε · T_1 · R(z, T_0))ⁿ`

converges (Neumann series) iff `‖ε · T_1 · R(z, T_0)‖ < 1`, giving

> `‖R(z, T_3)‖ ≤ ‖R(z, T_0)‖ / (1 − ε · ‖T_1‖ · ‖R(z, T_0)‖)`.

For this to be useful, need `T_0` with a **smaller** condition-number-and-explicit-resolvent and `T_1` small enough that the series converges.

## B.2 Natural splits

### B.2.1 Diagonal split: T_0 = diag(eigenvalues), T_1 = T_3 − diag

Take `T_0 := diag(1/2, 1/4, 1/8)` (a normal operator, explicit resolvent).

On γ: `‖R(z, T_0)‖ = max_λ 1/|z−λ| ≤ 8` (Approach A's spectral-radius bound, **exact** for the normal T_0).

`T_1 := T_3 − T_0` has entries:

> T_1 = ⎡ 7/8 − 1/2   −7/32       1/64     ⎤
>       ⎢ 1            0 − 1/4    0        ⎥
>       ⎣ 0            1          0 − 1/8  ⎦
>     = ⎡  3/8   −7/32    1/64 ⎤
>       ⎢  1    −1/4      0    ⎥
>       ⎣  0      1     −1/8   ⎦

This is **not small**: `‖T_1‖_F² = (3/8)² + (7/32)² + (1/64)² + 1 + (1/4)² + 1 + (1/8)² ≈ 0.141 + 0.048 + 0.0002 + 1 + 0.0625 + 1 + 0.0156 ≈ 2.27`. So `‖T_1‖_F ≈ 1.51` and `‖T_1‖_2 ≤ 1.51`.

Convergence condition: `‖T_1‖ · ‖R(z, T_0)‖ < 1`, i.e., `1.51 × 8 = 12.1`. **Series DIVERGES**.

The diagonal split fails because T_3 in companion form is not a small perturbation of its diagonalization — the off-diagonal `1`'s in the shift block are O(1), and on γ where `‖R(z, T_0)‖` is already 8, multiplication by an O(1) perturbation breaks Neumann.

### B.2.2 Shift-block split: T_0 = shift, T_1 = top row

Take `T_0 := [[0,0,0],[1,0,0],[0,1,0]]` (the nilpotent shift block). Then `T_0` has spectrum {0, 0, 0}; γ encloses none of these. `‖R(z, T_0)‖_op` on γ requires computing `(z·I − T_0)⁻¹` explicitly:

> `(z·I − T_0)⁻¹ = z⁻¹ · (I + (T_0/z) + (T_0/z)² + ...)`  if `|T_0/z| < 1`.

Since T_0 is nilpotent of index 3, `T_0³ = 0`, so the series terminates:

> `(z·I − T_0)⁻¹ = z⁻¹·I + z⁻²·T_0 + z⁻³·T_0²`.

On γ (|z| ≥ 1/2 − 1/8 = 3/8 at worst), `|z⁻¹| ≤ 8/3 ≈ 2.67`. So `‖R(z, T_0)‖ ≤ 2.67 + 2.67² + 2.67³ ≈ 27.8` (very loose, but T_0 is far from being "close to" T_3 — the perturbation is large).

`T_1 := T_3 − T_0 = ` `[[7/8, −7/32, 1/64], [0,0,0], [0,0,0]]`. `‖T_1‖_2 = ‖first row‖₂ = √((7/8)² + (7/32)² + (1/64)²) ≈ √(0.766 + 0.048 + 0.0002) ≈ 0.903`.

Convergence: `0.903 × 27.8 ≈ 25.1`. **Series DIVERGES**.

### B.2.3 Why neither split works

Both splits fail because **T_3's resolvent magnitude on γ is intrinsically large** (Approach A: M_3 in [8, 944]). Any base operator T_0 with computable resolvent will have its own R(z, T_0) of comparable magnitude on γ (because γ passes close to spec(T_3), and any "nearby" operator will have similar resolvent profile). The Neumann condition `‖T_1‖ · ‖R(z, T_0)‖ < 1` requires T_1 to be small relative to that intrinsic magnitude — but T_3's non-normality means there's no obvious small perturbation.

## B.3 Outcome

> **APPROACH_B_DIVERGES.**
>
> No natural split T_3 = T_0 + T_1 yields a convergent Neumann series for `R(z, T_3)` on γ. The diagonal split has perturbation `O(1)` vs base resolvent `O(8)`; the shift-block split has perturbation `O(1)` vs base resolvent `O(28)`. Both products `‖T_1‖ · ‖R(z, T_0)‖ ≫ 1`.

Approach B does not refine the bound from Approach A.

## B.4 Honest scope

(1) Other splits (e.g., perturbation around a 2×2 block) might converge, but the natural choices fail.

(2) Approach B is a refinement tool; it requires Approach A's bound to be loose by a specific factor that a perturbation argument can recover. Approach A's bound is already inside R77.2's stated range; no obvious slack for Approach B to exploit.

(3) As with Approach A, this calc is moot under R77.3 falsification — even if Approach B converged, it would refine a number computed for a matrix that doesn't describe ε_n.
