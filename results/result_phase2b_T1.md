# Result — THEOREM (Real-T1) PROVEN + gate GREEN: the real q=3 pair-operator eigenvalues are the twisted autocorrelations of the halving weights. The kinematic half of the q=3 boundary is closed.

**Date:** 2026-07-16. Wilson's closed-form theorem for the real q=3 pair operator, with the pre-registered agent gate (18-member family at L=3, exact ≤1e-12) **PASSED**. This is the program's **second proven result — and the first on the real operator.** Probe `probes/probe_phase2b_T1.py`, log `logs/probe_phase2b_T1_log.txt`.

**Headline: the eigenvalues of the real Collatz pair operator ARE the twisted autocorrelations of the halving weights; the eigenvectors are their correlation profiles, supported entirely on zero carry. No dressing, no resolvent, no invariant hunt. Gate: L=2 6/6 (worst 7.08e-16, convention locked), L=3 18/18 (worst 3.77e-15) — all ≤ 1e-12.**

## THEOREM (Real-T1; Wilson). Kinematic half of the q=3 boundary.
For the pair operator at q=3, any L, with `D = 2·3^{L−1}` and folded weights `w_δ ∝ 2^{−δ}`: for each `k ∈ ℤ/D`, the function
```
ℓ_k(a, b, γ) = ω^{−e_a} · R_k(e_ρ)/R_k(0) · [γ = 0],
    ω = e^{2πik/D},   R_k(e) = Σ_δ w_δ w_{δ+e} ω^δ,
    e_a = log₂(a),   e_ρ = log₂(b·a⁻¹)
```
is an **exact left eigenvector** with eigenvalue `c_k = R_k(0) = Σ_δ w_δ² ω^δ`.

**Proof (three one-line steps):**
1. **Protection** — γ′ = 0 forces γ = T = 0 in representative arithmetic (nonnegative integers can't cancel), so the eigen-equation holds trivially off the zero-carry sector.
2. **Funnel** — on {γ=0}, T = 0 exactly forces a′ = b′, so the visible dynamics per twist sector is rank one: (Mᵀf)(e) = R_k(e)·f(0).
3. **Rank-one spectrum** — unique nonzero eigenvalue R_k(0), eigenvector = the autocorrelation profile. ∎

The exact eigenvalues are the twisted autocorrelations `R_k(0)` of the halving weights; the eigenvectors are their correlation profiles `R_k(e)/R_k(0)`, carried on {γ=0} with the twist ω^{−e_a}.

## Gate — 18-member family at L=3 (pre-registered exact ≤ 1e-12). ✅ PASS
Method: build M = `build_M_gen(3, L, 2, [λ^δ])`, construct ℓ_k from the formula, measure `‖Mᵀℓ_k − c_k ℓ_k‖∞ / ‖ℓ_k‖∞`. Direct sparse mat-vec (no ARPACK, per the instrument law).

| level | members | worst relative residual | verdict |
|---|---|---|---|
| L=2 (D=6) | 6/6 | 7.08e-16 | convention locked (Wilson's 6/6 reproduced; s=+1, first try) |
| **L=3 (D=18)** | **18/18** | **3.77e-15** | **ALL PASS (≤1e-12)** |

`R_k(0) = c_k` are exactly the circulant-family eigenvalues (already H_CIRC-confirmed 18/18); this gate confirms the **closed-form eigenVECTORS** at every k, on the full dressed support (all e_ρ, including the off-diagonal and R1's split-class states). **Claude-verified.**

## Three inversions (banked)
1. **The obstruction IS the engine.** The floor-arithmetic that killed S, J, and lumpability — the thrice-measured **0.258 obstruction** — is exactly the **protection lemma** (γ′=0 ⟹ γ=T=0). We spent three probes fighting the mechanism that proves the theorem.
2. **The STOP discipline caused the discovery.** The planned invariant hunt over "21 classes" is **moot — superseded, not failed.** Had R1 fished for a gauge yielding 21 clean classes, we'd have hunted a digraph invariant that doesn't exist as posed, while the truth (six support classes, rank one, orthogonal to that framing) sat unseen. **The refusal to decorate caused the discovery.**
3. **D1's toy architecture was the real structure — co-variant, not invariant.** The "invariant ray + solvable block" of D1 is here the **zero-carry sector (the ray) + a rank-one block**; the covariance is the ω^{−e_a} twist.

## Retractions / corrections
- **Walk-back #15 (Wilson):** "126 = 6·21 ⟹ 21 clean classes" was divisibility numerology (inferring gauge-invariance from a factorization coincidence). Retracted. R1's STOP stands; the closed form supersedes the class reduction entirely.
- **L-B evidence strengthened.** L-B (LALB) verified k=0 gauge-invariance — but for k=0 the twist ω^{−e_a} is trivial (ω=1), so that check was weak/near-vacuous and never touched the dressed off-diagonal support for k≠0. This gate verifies the **full closed form for all 18 k on all states** — strictly stronger. (State-level reachability splitting, R1, does not refute the function-level reduction: the observable dynamics descends per twist sector, which is what the eigen-equation needs.)

## Board after this turn
- **Phase 2's kinematic half is DONE** (theorem proven in sketch, L-independent arithmetic; gate green at L=2, L=3) — pending Wilson's write-up. **Bankable as the program's second proven result, first on the real operator.**
- **The one remaining object of the entrance exam: the DYNAMICAL PARTNER** — the braiding mode, the rate law (2.9e-3, 1.0e-4, ?), the coalescence with `c₀ = R(0)` as L→∞. D3 now has a **proven family to measure the partner against**: the partner is *the one thing in the spectrum near 1/3 that is not an autocorrelation.* R2 (L=4 partner) remains the standing target for that (deferred: local wall + Lambda option).

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, the H-gates, the J-refutation, L-A. No `r_q` value changes. R1's STOP is preserved (now superseded, not pending a gauge).

_Reporting discipline: the theorem is Wilson's (banked with its sketch proof); the gate is Claude-verified — pre-registered exact ≤1e-12, locked on L=2's 6/6 before the L=3 read, 18/18 pass at worst 3.77e-15 via direct sparse mat-vec. The invariant-hunt cancellation and walk-back #15 are recorded as supersession, not failure. L-B's prior evidence is honestly re-scoped as weaker than the full closed form now verified._
