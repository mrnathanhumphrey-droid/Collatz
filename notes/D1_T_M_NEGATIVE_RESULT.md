# D1 (full T_M on V_n^M / bilinear V_n ⊗ V_n*): NEGATIVE RESULT

**Date:** 2026-05-15. Follow-up to K_k structural lemma + U_n / Phi_omega probe. Tests whether the bilinear pair-correlation T_M operator has discrete eigenvalues matching T_lead's 43/45 (Q-rigorous within-level rate) and/or empirical period-9.2 CC pair structure.

## Headline

> **The bilinear T_M operator (in BOTH attempted forms — V_n^M truncation and full V_n ⊗ V_n* tensor) has spectrum max |λ| ≈ 0.345 at n=5, far below T_lead's 43/45 ≈ 0.956 or empirical ρ ≈ 0.984.** No CC pair structure emerges matching period-9.2.
>
> **T_lead's 43/45 is a CLASS-RESOLVED COHERENT-SUMMATION phenomenon**, NOT an eigenvalue of any natural bilinear pair-correlation operator at the full tensor or V_M level. The (P_+, P_-) projection achieves phase-coherent amplification of small individual contributions, similar to how Plancherel mass coherently combines from individual ξ-modes.

## Two attempts, both confirming structural mismatch

### Attempt 1: Truncated T_M on V_n^M (η-parameterized moments)

Built T_M_trunc as a linear operator on V_n^M (dim N_n = 2·3^{n-1}) via Tao recursion with the δ ≡ 0 (mod 3^{n+1}) truncation (the "drop phase-twist contributions" projection).

Result: spectrum at n=2, 3, 4 has a single non-zero eigenvalue at ≈ 1/3, all others zero. (P_+, P_-) projection gives 0.185, 0.062, 0.021 at n=2, 3, 4 — **NOT 43/45**.

The truncation is too aggressive: it drops the η' ≡ 1 (mod 3) but η' ≠ 1 contributions that T_lead's T_diag uses (where the phase character factors through level n via 3 | δ, not necessarily 3^{n+1} | δ).

T_lead lives on a different structural object: the cross-frequency moment space V_M with class-resolved (P_+, P_-) parameterization, not V_n^M.

### Attempt 2: Bilinear tensor U_n ⊗ conj(U_n) on V_n ⊗ V_n*

Built the natural bilinear pair-correlation operator at the tensor product level. The pair-correlation Tao recursion lifts to U_n ⊗ conj(U_n): V_n ⊗ V_n* → V_{n+1} ⊗ V_{n+1}*. Using twisted-fiber projection on both factors gives endomorphisms on V_n ⊗ V_n*.

Tensor product spectra (top |z| across n=2..5):

| Operator | n=2 | n=3 | n=4 | n=5 | Limit |
|---|---|---|---|---|---|
| Phi_ω ⊗ conj(Phi_ω) | 0.277 | 0.334 | 0.101 | 0.102 | ~0.102 = 0.319² |
| Phi_ω² ⊗ conj(Phi_ω²) | 0.277 | 0.334 | 0.346 | 0.345 | ~0.345 = 0.587² |
| Mixed Phi_ω ⊗ conj(Phi_ω²) | 0.155 | 0.182 | 0.187 | 0.187 | ~0.187 = 0.319·0.587 |

**All tensor moduli converge to PRODUCTS of Phi_omega eigenvalues** (which were 0.319 and 0.587 from the previous probe). No eigenvalue near 43/45 = 0.956 anywhere in the spectrum.

Closest to 43/45 across all three tensors: 0.345 (Phi_ω² ⊗ conj(Phi_ω²) at n=5). 64% gap. Not a candidate.

**Mixed-tensor eigenvalues have non-trivial angles**: at n=5, the top |z| ≈ 0.187 eigenvalues have args ∈ {-0.46, -0.51, -0.53, +2.63, ...} rad, continuously distributed in arcs, NOT concentrated at θ = 2π/9.2 ≈ 0.683. **Continuous-spectrum-on-arc structure, no discrete CC pair.**

## What this means structurally

T_lead's 43/45 emerges as **1 − Σ_g W_+(g) = 1 − 2/45** through a specific class-resolved coherent summation:
- T_diag = (1/5)·[[1,1],[4,4]] has eigenvalue 1 on (1, 4) (R64.B class-mass + Plancherel)
- Off_lin contributes −2/45 on (1, 4) (cross-frequency closure at g ∈ {2, 4, 6, ...})
- T_lead = T_diag + Off_lin has eigenvalue 43/45 on (1, 4)

This 43/45 is NOT visible at the bilinear pair-correlation tensor level (max |eig| ≈ 0.35) because that level sees individual pair contributions of order Phi_omega's eigenvalues (≈ 0.32, 0.59), and the **phase-coherent sum over the (P_+, P_-) projection** is what amplifies these into 43/45.

It's analogous to Plancherel mass: individual Fourier coefficient magnitudes |μ̂_n(ξ)| are small (≈ 1/3^{n/2}), but the SUM Σ |μ̂_n(ξ)|² = S_n approaches 7/15. The coherent sum is the rate carrier; individual modes aren't.

## Why this is consistent with R77.6 + T_V_DISPOSITION

R77.6 found E(z) has a branch cut at z=2 (continuous spectrum endpoint). T_V_DISPOSITION found V_M doesn't close under iteration (phase offsets + odd-G shifts). Both indicate: **finite-truncation discrete operators DON'T carry the c=7/45 asymptotic rate as a single eigenvalue.**

My probes (K_k, U_n, Phi_omega, tensor T_M) all confirm:
- K_k: spectrum {1, 0, ..., 0}, no rate info
- U_n: pair-structured σ in 0.77..1.29 (non-trivial but not c=7/45)
- Phi_omega: continuous-on-circle at radii 0.319, 0.587
- Tensor T_M: max |eig| ≈ 0.345, no 43/45 anywhere

**Conclusion: the c=7/45 rate (whether 43/45, 0.984, or other) is fundamentally a CLASS-RESOLVED COHERENT-SUMMATION phenomenon at the (P_+, P_-) reduction, NOT an eigenvalue of any natural finite-truncation operator on V_n / V_n^M / V_n ⊗ V_n*.**

## Implications for the period-9.2 CC pair

The period-9.2 oscillation in the empirical sign pattern (+ + − − − − − − − + + + +) at n=2..13 is NOT found as a discrete CC pair in any of the inter-level operators tested. Either:

**Option I.** The period-9 is a finite-n transient that doesn't persist to the asymptotic (consistent with PADE_NUMERICAL_DISPOSITION's "transient ratios at n=10..13").

**Option II.** The period-9 is a continuous-spectrum resonance at a Mellin saddle point — not visible at finite n in operator eigenvalues, but real in the analytic continuation of E(z).

**Option III.** The period-9 is a phase-coherent class-resolved phenomenon at a NON-(P_+, P_-) projection — analogous to how 43/45 lives on (P_+, P_-). Identifying this projection would surface the CC pair.

Option III is the most actionable. It would require constructing a class-resolved projection of T_M that produces a 2x2 (or larger) operator whose spectrum DOES include a CC pair. This is similar to how T_lead emerges from (P_+, P_-) projection; the period-9 CC pair would emerge from a different specific projection.

The natural candidates for "different class-resolved projection":
- (P_++, P_+-, P_-+, P_--) 4-dim class-pair-resolved
- Higher-order class moments (P_+++ etc.) at 8-dim and beyond
- Phase-resolved projections (P_+,φ, P_-,φ for specific phase parameters φ)

This is the substantive next probe if pursued. Cost estimate: 1-3 sessions.

## Files

- [`T_M_truncated_spectrum.py`](T_M_truncated_spectrum.py) — Attempt 1 (V_n^M truncation)
- [`T_M_tensor_spectrum.py`](T_M_tensor_spectrum.py) — Attempt 2 (tensor V_n ⊗ V_n*)
- [`experiments_output/T_M_truncated_spectrum.json`](experiments_output/T_M_truncated_spectrum.json)
- [`experiments_output/T_M_tensor_spectrum.json`](experiments_output/T_M_tensor_spectrum.json)
- [`D1_T_M_NEGATIVE_RESULT.md`](D1_T_M_NEGATIVE_RESULT.md) — this writeup
- Companion: [`T_V_DISPOSITION.md`](T_V_DISPOSITION.md) (V_M doesn't close — algebraic version)
- Companion: [`INTERLEVEL_U_PROBE_RESULT.md`](INTERLEVEL_U_PROBE_RESULT.md) (Phi_omega continuous-on-circle)
- Companion: [`T_LEAD_CORRECTED_DISPOSITION.md`](T_LEAD_CORRECTED_DISPOSITION.md) (43/45 via class-resolved 2x2)
