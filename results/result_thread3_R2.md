# Probe R2 — the amplitude (the 7/10 gate; thread 3 session two)

**Date:** 2026-07-19  CPU, exact. Probe `probes/probe_thread3_R2.py`. **A∞ = 7/10 holds algebraically; the
crown's stated closed-form route (the ℓ₀ overlap) FALLS, for a provable reason — reported straight.**

Object algebra (pinned): A_k := 3^k·P(agree mod 3^k) = 3^k·Σ_{m≥k} g_m, g_m = a_m − a_{m−1}/3, a_m = ‖π_m‖².
Recursion **S_k = A_k − A_{k+1}/3 = 3^k g_k = c_k** (R1's welded shell). Since c_k → 7/15, **A∞ = (3/2)(7/15) = 7/10**.

## R2-A — accumulation object + consistency: **CONFIRMED**
- **Exact consistency (the weld):** S₁ = 3g₁ = **2/3**, S₂ = 9g₂ = **10/21** (exact rationals, `==` True). The
  recursion S_k = A_k − A_{k+1}/3 reproduces the R1 welded shell sequence exactly.
- **Algebraic crown value:** A∞ = (3/2)·S∞ = (3/2)(7/15) = **7/10 exact**.
- **Measured** A_k (frozen L=3) → ~0.66–0.69, i.e. **(3/2)×(L=3 shell level ≈0.44)** — the same −2.5% shell
  truncation as R1 (A_10 = 0.6688 ≈ (3/2)·0.446). Approaches 7/10 as L→∞ at the R73 rate. Object confirmed.

## R2-B — the crown's ℓ₀ ansatz: **DOES NOT evaluate to 7/10 (and the reason is a one-line theorem)**
Wilson's pre-registered closed form was ⟨ℓ₀|v_indep⟩ with ℓ₀ = R₀(e_ρ)/R₀(0) (Real-T1's zero-carry c₀ left
eigenvector) and v_indep = w⊗w. This reduces exactly to **Σ_e R₀(e)²/R₀(0)** and evaluates to:

| L | ⟨ℓ₀\|v_indep⟩ (exact) | float | vs 7/10 |
|---|---|---|---|
| 2 | **2531/4095** | 0.61807 | −0.082 |
| 3 | 12726877601/22906492245 | 0.55560 | −0.144 |

**Not 7/10, and trending AWAY (0.618 → 0.556)** — reported as the exact rational it *does* evaluate to, per the
spec. This is not finite-D noise (R₀(0)=Σw² is already 0.33334 at L=3); the functional is structurally wrong.

**Why (airtight, one line):** the c₀ **simple** mode is *invisible* to the shell functional. Its contribution to
g_m = a_m − a_{m−1}/3 is A_{c₀}[(1/3)^m − (1/3)·(1/3)^{m−1}] = A_{c₀}[(1/3)^m − (1/3)^m] = **0**, exactly, ∀m. So
A∞ receives **nothing** from the isolated c₀ eigenvector — the amplitude lives entirely in the **partner-Jordan
coupling** at the degenerate 1/3 (c₀ and partner coincide there — the crossing): **A∞ = (3/2)·β**, β = the
k·(1/3)^k coefficient of a_k = the Jordan slope = S∞ = 7/15.

**Correction for the pen (walk-back on the crown's route, not the value).** 7/10 is *not* the naive kinematic
c₀-overlap; it is the **generalized (Jordan) eigenvector coupling** of the degenerate 1/3 eigenspace. The closed
form must use the generalized left eigenvector g₀ (the tower-coupled Jordan partner of c₀), not ℓ₀ alone — since
the crossing makes c₀ defective, "overlap with the c₀ mode" must mean the whole 1/3 Jordan block. With that
substitution the amplitude is (3/2)·β = 7/10 whenever β = 7/15 (the flat-level theorem). **Caveat on
independence:** routed through β = S∞, the 7/10 is then *algebraically equivalent* to the flat-level theorem, not
a fully independent third derivation — unless the Jordan slope β is derived from the Jordan structure without
invoking the renewal limit. The isolated-ℓ₀ route that *would* have been independent is the one that fails here.

## R2-C — q=5,7 contrast: **A_k → 0 (agreement amplitude vanishes)**
| q (L) | A_1..A_10 | behavior |
|---|---|---|
| 5 (L=2) | 0.149, 0.048, 0.027, 0.019, 0.015, 0.010, 0.004, 0.0002, −0.001, −0.001 | → 0 |
| 7 (L=2) | 0.305, 0.058, 0.033, 0.016, −0.003, −0.007, 0.002, … | → 0 |
A_k decays toward 0 — the spectral gap seen as vanishing agreement amplitude (no persistent c₀/partner degeneracy
to sustain it). The precise geometric ratio (3/q) is contaminated at L=2 (A_k hits the O(1e-3) band-noise and
sign-flips by k≈8, same L-truncation as R1's c_k) — vanishing is the clean signal, the ratio is not.

## Status
**R2-A confirmed** (exact welds 2/3, 10/21; A∞ = (3/2)·S∞ = 7/10 algebraic). **R2-B: the crown's isolated-ℓ₀
closed form evaluates to 2531/4095 (L=2), NOT 7/10 — the c₀ simple mode is provably invisible to the shell
functional, so 7/10 is the degenerate-1/3 Jordan coupling, needing the generalized eigenvector.** R2-C: q=5,7
agreement amplitude → 0. The bridge's value (7/10 → 7/15 → 7/45) stands; the crown's *route* (kinematic c₀
overlap) is corrected to the Jordan/tower coupling, with a flagged caveat on its independence. No fitting;
exact rationals reported.
