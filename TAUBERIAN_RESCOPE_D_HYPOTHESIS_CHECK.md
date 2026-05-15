# TAUBERIAN_RESCOPE_D_HYPOTHESIS_CHECK (Newman-Zagier Analytic Theorem × inputs)

**Date:** 2026-05-13.

## h × I matrix

Try construction T1: f(t) = ε_⌊e^t⌋ · e^{-t} (some normalization on log timescale).

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1: f bounded on [0, ∞) | UNVERIFIABLE for full sequence; empirically |ε_k|·2^k ≤ 0.4 ≈ bounded, so f(t) = ε_⌊e^t⌋ would be bounded if ε_n bounded. PLAUSIBLY SATISFIED. | | | |
| h_2: f locally integrable | SATISFIED for any reasonable construction. | | | |
| h_3: g(z) := ∫_0^∞ f(t) e^{-zt} dt exists for Re z > 0 | Reduces to Σ ε_n e^{−z log n} = Σ ε_n n^{-z} (a Dirichlet series), convergent for Re z > σ_a some abscissa. PLAUSIBLY SATISFIED if ε_n bounded. | | | |
| h_4: g extends holomorphically to Re z ≥ 0 | UNVERIFIABLE — Mode H circular. **Same load-bearing condition as C: g's holomorphic extension to Re z ≥ 0 is what we'd derive from a closure result, not assume.** | | | |

Output for SELECTED case would be: ∫_0^∞ f(t) dt converges and equals g(0). This is even *weaker* than Wiener-Ikehara — only convergence of an integral, no asymptotic.

For the c=7/45 closure, we'd need to identify g(0) = 7/45 (or a rational related to 7/45). Without h_4, we cannot evaluate g(0).

**Theorem D disposition: BLOCKER (h_4 Mode H circular).**

---

## Observation: D is strictly weaker than C

Newman-Zagier Analytic Theorem requires h_4 (g extends holomorphically to Re z ≥ 0) — the closed right half-plane. Wiener-Ikehara C 4.2 requires the analogous condition on Re z = 1 (the abscissa of convergence shifted). Both conditions are Mode H circular for Syracuse.

D is more general (it doesn't assume monotone nondecreasing), but the analytic-continuation condition is identical in *flavor*: the Dirichlet/Laplace transform must extend holomorphically past the abscissa of convergence — a property that is essentially the closure target itself.

**Aggregate disposition for D: BLOCKER (Mode H circular).**
