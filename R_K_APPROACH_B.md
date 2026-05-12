# R_K_APPROACH_B — perturbation series / Neumann decomposition for R_k

**Date:** 2026-05-11. Phase 2B of R_K probe. Considers whether a Neumann-series decomposition R_k = R_0 + δR around a known-norm baseline can bound an effective resolvent on R_k.

## 1. The structural obstruction

Approach A established (R_K_APPROACH_A.md §6) that the candidate finite-dimensional Φ_k : W_{k−1} → W_k built from (lift-then-residual) **does not actually transport R_{k−1} → R_k** — the regression coefficient c_k = 0 exactly (R77.5 §3, structural from marginal consistency). So R_k is **not in the image of any natural Φ_k from a previous level**; it's a fresh perpendicular contribution at each level.

Approach B asks: can we still find a baseline "R_0" with known norm, write R_k = R_0 + δR with ‖δR‖ small, and run Neumann to control a resolvent?

## 2. Candidates for R_0

### 2.1 Cardinality-uniform R_0

The leading-order behavior of R_k is the cardinality-uniform contraction ‖R_k‖² = c · 3^{-k} with c ≈ 0.155.

Natural "R_0" candidates:

- **Constant-magnitude random vector in W_k**, normalized to ‖R_0‖² = c · 3^{-k}.
- **R_k itself at level k − 1** lifted to level k (T(R_{k−1})): but ‖T(R_{k−1})‖² = ‖R_{k−1}‖²/3 = (c · 3^{-(k-1)}) / 3 = c · 3^{-k}, the right norm — and ⟨R_k, T(R_{k−1})⟩ = 0 exactly. So δR = R_k − T(R_{k−1}) has ‖δR‖² = ‖R_k‖² + ‖T(R_{k−1})‖² = 2c · 3^{-k} = 2 ‖R_0‖² > ‖R_0‖². Neumann fails immediately: ‖δR‖/‖R_0‖ > √2 > 1.
- **The c-uniform vector itself** (the 7/45-rate prediction with no structural fluctuation): R_0 = c · 3^{-k/2} · u for some unit u. But R_k − R_0 captures all the variation, and that variation has ‖R_k − R_0‖² of the same order ‖R_k‖² ≈ c · 3^{-k}. Again ratio O(1), Neumann fails.

### 2.2 Why all decompositions of R_k blow up the perturbation

R_k has **no leading low-rank structure that survives across k**. Two observations:

(a) **R_k's effective rank grows with dim(W_k) = 4·3^{k−1}.** The per-coordinate squared mass ‖R_k‖²/dim(W_k) ≈ 0.116 · 9^{-k} stays bounded; mass is **spread across all coordinates of W_k**, not concentrated in a low-dim subspace.

(b) **Cross-level correlation is zero structurally:** c_k = 0 over Q (R77.5 §3). So no rank-1 update of T(R_{k−1}) captures R_k; the closest rank-1 approximation in T(W_{k−1}) direction is the zero vector.

Together: R_k is **structurally close to a maximally-spread (incoherent) vector in W_k** with the L² norm carrying all the information. There is no "R_0 with small perturbation" structure to exploit.

## 3. Operator-level Neumann (for the would-be Φ_k)

If we ignore the c_k = 0 obstruction and ask "can we write Φ_k = Φ_0 + δΦ with ‖δΦ‖ < something / ‖(zI − Φ_0)^{-1}‖":

From Approach A §4.1 the natural prediction is σ_1(Φ_k) → 1 as k grows. Two splits:

### 3.1 Split via top singular value

Φ_0 := σ_1 · u_1 v_1^T (rank-1, leading mode).
δΦ := Φ_k − Φ_0 (bulk modes).

- ‖Φ_0‖ = σ_1 ≈ 0.75..0.91 (k=2..5).
- ‖δΦ‖ = σ_2 ≈ 0.62 (bulk).

Contour γ around z = 1/2 in σ²-space:

- distance from 1/2 to σ_1² ≈ 0.56..0.83 is ≈ 0.06..0.33 (small, especially as k grows).
- ‖(zI − Φ_0^* Φ_0)^{-1}‖ = max(1/|z − σ_i²|, 1/|z|) on γ = 1/(0.06..0.33) ≈ 3..16 at the closest approach (depends on contour radius).
- Neumann condition: ‖(zI − Φ_0)^{-1}‖ · ‖δΦ‖ < 1, requires ‖δΦ‖ < 0.06..0.33, but ‖δΦ‖ ≈ 0.62 ≫ that. **Neumann fails.**

### 3.2 Split via near-zero bulk

Φ_0 := bulk modes (σ_2, σ_3, ..., σ_min), Φ_1 := σ_1 · u_1 v_1^T.

- ‖Φ_0‖ = σ_2 ≈ 0.62.
- ‖Φ_1‖ = σ_1 → 1.
- Same problem: bulk's "spectrum" (the σ_i² ≈ 0.38..0.44) is near 1/2 = 0.5, but contour-resolvent norm relative to σ_2 needs the contour to separate from 0.38, which gives ‖(zI − Φ_0)^{-1}‖ ≈ 1/0.12 ≈ 8, and ‖Φ_1‖ → 1 dominates. **Neumann fails.**

## 4. Why this echoes the M_3 probe's Approach B

The M_3 probe (`M3_APPROACH_B.md`) tried the same Neumann decomposition on T_3 (companion form) and found both natural splits diverge. The structural reason there: non-normality of T_3 forces κ(V) to be the limiting factor, and the perturbation δT_3 wasn't small relative to (zI − T_0)^{-1}'s norm.

Here the reason is different but analogous: **σ_1(Φ_k) → 1 means the leading mode dominates and the bulk near σ² ≈ 1/2 is not isolatable**. No clean rank-1 + small-perturbation decomposition exists.

## 5. Outcome

**APPROACH_B_FAILS.**

Both natural splits of Φ_k (top-mode-isolation and bulk-isolation) give Neumann ratios well above 1. **No effective M_3' bound** via perturbation.

This is consistent with: **Φ_k is not the right operator to bound** (because c_k = 0 means it doesn't carry R_{k−1} → R_k anyway). Approach B's failure is downstream of Approach A's structural finding, not an independent obstruction.

## 6. What would an "R_0" baseline have to look like to succeed?

To get a clean Neumann, we'd need an operator Φ_0 such that:

- ‖Φ_0‖ < 1/2 (or some explicit value < 1).
- Φ_0 has a clean resolvent in a contour separated from σ²-target.
- δΦ := Φ_k − Φ_0 has ‖δΦ‖ · ‖(zI − Φ_0)^{-1}‖ < 1 uniformly in k.

The natural candidate "K_k as a contraction with eigenvalue 1 isolated" doesn't work because K_k|_{W_{k−1}} (the restriction to mean-zero-fiber subspace) has ‖K_k|_{W_{k−1}}‖ near 0 already, not near 1/2 — it's a **strongly mixing** operator with no spectrum near 1/2 (R77.4 erratum).

There is **no project-characterized operator** whose spectrum has the rate-1/2 feature with a resolvent contour around 1/2 isolating it. This is the same obstruction the M_3 probe hit (M3_DISPOSITION.md).

## 7. Files

- `R_K_APPROACH_B.md` (this file) — perturbation outcome
- `M3_APPROACH_B.md` — parallel outcome for T_3 (same FAILS verdict, different reason)

Phase 2C (resolvent-norm numerical) follows but is fundamentally moot for the same reason: there is no canonical operator to compute the resolvent of, because the candidate Φ_k doesn't actually carry the R_{k−1} → R_k dynamics.
