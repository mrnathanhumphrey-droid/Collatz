# Probe R17 — the slow mode — **A convention (with a coordinate caveat), D PASS, B/C quasi-stationary**

**Date:** 2026-07-21  Reuses R7/R10. Probe `probes/probe_slowmode_R17.py`. Studies the transport symbol D and the
deviation field δ_r in Fourier, following R16's transport gate.

## R17-A — D CONVENTION (forced): the closed form is the *linearized* slow-mode symbol
The proposed **|D(ξ)|² = 1/(5 − 4cos(πξ/3^r))** is the symbol of the idealized slow-mode shift D_ideal(ξ) =
Σ_{v≥1} 2^{−v} e(−ξv/(2·3^r)) (the dlog shift −v·dlog₄(2) = −v/2, geometric series — exact). For this D:
- **(i) conjugation |D(ξ)| = |D(−ξ)| — exact** (dev ≤3e-16, all r).
- **(ii) circle-average of |D_ideal|² = 1/3 exactly** (measured 0.333333 to 6 digits, r=3…7) = Σ_v p_v². ✓
- **(iii) holds by construction** (it *is* the closed form).

**Coordinate caveat (honest, flagged for the pen):** the *exact* transport (R16-A) is **not diagonal** — in
additive-Y Fourier it is μ̂_r(ξ) = E_v[e(ξ·2^{−v}/3^r)·μ̂_{r−1}(ξ·2^{−v})], where the frequency **scales**
ξ → ξ·2^{−v}, not a pure shift. The directly-computed v-averaged additive phase D_add(ξ) = E_v[e(ξ·2^{−v}/3^r)]
is a *different* object: it satisfies (i) exactly, gives (ii) = 1/3 + O(2^{−ord}) (65/189 at r=2, machine-1/3 by
r=4), and **does not** match (iii) (deviations 0.5–0.9). So "D" with the proposed closed form is the **leading /
linearized slow-mode symbol**, correct as such; it is not the exact transport symbol (there is no clean diagonal
one). #37 not incurred — (i)/(ii)/(iii) all hold for the intended D_ideal; the exact-transport non-diagonality is
noted as scope, not a failure.

## R17-D — FIXED-POINT FORM (forced): **PASS**
T is a self-map — 1 + 3·2^{−v}·X ∈ 1+3ℤ/3^{r+1} for all v (X≡1 mod 3 ⟹ 1+3·unit ≡ 1 mod 3), r=2…5 — and the
achieved-resolution invariance holds: **μ_{r+1} folded mod 3^r == μ_r** exactly (Tao 1.22 consistency), r=2…5. The
fixed-space reframe is certified as computation.

## R17-B — QUASI-STATIONARY SELF-CONSISTENCY (measurement, NO fit): δ contracts at the mean rate 1/3
δ-weighted mean contraction ⟨|D|²⟩_{δ_r} := Σ_ξ |δ_r(ξ)|²|D(ξ)|² / Σ_ξ|δ_r(ξ)|² (δ_r(ξ)=μ̂_r(ξ), ξ primitive):

| r | ⟨\|D_ideal\|²⟩_{δ_r} | ⟨\|D_add\|²⟩_{δ_r} | plain avg \|D_ideal\|² |
|---|---|---|---|
| 3 | 0.305742 | 0.372761 | 0.333332 |
| 4 | 0.334466 | 0.415925 | 0.333333 |
| 5 | 0.332214 | 0.426798 | 0.333333 |
| 6 | 0.334300 | 0.428232 | 0.333333 |
| 7 | 0.333887 | 0.429149 | 0.333333 |

**⟨|D_ideal|²⟩_{δ_r} ≈ 1/3, flat in r (r≥4), tracking the unweighted circle-average exactly.** The deviation field
contracts at the **mean rate 1/3** — it is *not* concentrated in a slow region (that would give ⟨|D|²⟩ > 1/3) nor
anti-concentrated (< 1/3); it sits at the average. This is the quasi-stationary balance: δ's magnitude contracts by
~1/3 per step but is replenished by the DC source, holding ‖δ_r‖² ≈ S_r → 7/15 at steady state. (The 1/3 is the L²
contraction; the smallness/oscillation of the *signed* bulk b_r is a phase effect on top of this magnitude balance,
consistent with R15/R16.) The D_add column (drifting to ~0.429) is the different additive-phase object, reported for
completeness; the slow-mode symbol is D_ideal. Per-r statistic, no cross-r rate extracted.

## R17-C — LOCALIZATION SCALE (measurement, NO fit): the angular width HOLDS
|δ_r|²-weighted moments of x = dist(ξ/3^r, 0):

| r | ⟨x⟩ | √Var | frac ‖δ‖² in x<1/12 |
|---|---|---|---|
| 4 | 0.25305 | 0.15144 | 0.1500 |
| 5 | 0.24751 | 0.13910 | 0.1616 |
| 6 | 0.24915 | 0.14340 | 0.1318 |
| 7 | 0.24944 | 0.14435 | 0.1702 |

**The mode's angular width holds (√Var ≈ 0.14, stable; ⟨x⟩ ≈ 0.25, stable) — quasi-stationary, neither shrinking
(being squeezed into a slow region) nor spreading (leaking).** The deviation field is **broadly spread** (⟨x⟩ ≈ 1/4,
the value for a near-uniform spread over the half-circle), not localized near x=0 — consistent with R14's ruling
that no non-uniform limiting profile exists. No mechanism claimed; the fraction in x<1/12 (~0.13–0.17) is
depth-noisy.

## R17-E — outstanding (one line each)
- **Symbol identity χ(4)·w(χ) = 1/(4 − χ(2)^{−2})** — verified numerically (r=3, k=1,2,4): χ(4)w(χ) =
  χ(4)/(4χ(4)−1) = 1/(4 − χ(4)^{−1}) and χ(2)^{−2} = χ(4)^{−1} ✓ (the weight is the halving-character symbol).
- **R85 rung-1 (n=8):** one Bluestein/support-pruned probe — DEFER (unchanged from R13-E/R16-E).
- **Branch weights / λ-sweep:** not specced here; carried as owed, no run.

## Status
**R17: A convention** (the closed form 1/(5−4cos(πξ/3^r)) is the linearized slow-mode symbol D_ideal — (i)
conjugation exact, (ii) circle-average = 1/3 exact, (iii) by construction; the *exact* transport is non-diagonal,
so the direct additive phase D_add differs — flagged as scope, #37 not incurred), **D PASS** (self-map +
achieved-resolution invariance), **B quasi-stationary** (⟨|D|²⟩_{δ_r} ≈ 1/3 flat — δ contracts at the mean rate,
balanced by the source; a QSD, not a localized slow eigenmode), **C width holds** (angular width ≈0.14 stable,
broadly spread not localized — consistent with R14 no-non-uniform-ψ). The deviation field is a broad quasi-stationary
field contracting at 1/3 in L²; the theorem's difficulty is the *signed* phase residual (bulk b_r), not the
magnitude. **Still owed (pen):** the uniform tower-contraction estimate on the signed field (the R16/R5 shared
crux). No fitting; exact self-map/identity gates, labeled numeric symbol/QSD/localization; the D coordinate
subtlety reported honestly.
