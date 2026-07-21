# Probe R9 — the collision identity gate + γ-tables — **FULL PASS (A gate + B/C/D)**

**Date:** 2026-07-21  Exact rationals. Probe `probes/probe_gamma_R9.py` (reuses R7 build_mu/mu1/cram). Gates
Wilson's session: the exact-at-every-level restatement

> **S_K = 2·Σ_{m≥1} 4^{−m}·γ_{K−1}(τ_m)**,  τ_m := (4^{−m}−1)/3 ∈ ℤ₃ (v₃(τ_m)=v₃(m)),
> **γ_n(τ) := 3ⁿ·Pr_{a,a′~μ_n iid}[(a−a′)+τ(1+3a) ≡ 0 mod 3ⁿ]**,  γ_n(0)=3ⁿ‖μ_n‖²=X_n.

γ computed **directly from the μ tables** (via the unique-collision partner a′=a+τ(1+3a) mod 3ⁿ, O(support)/τ),
**not routed through the C-tables** — the identity is tested, not assumed.

## R9-A — COLLISION IDENTITY GATE: **PASS**
S_K from the γ-engine = **frozen S_K exactly, K = 2, 3, 4, 5, 6** (byte-equal rationals). Hand anchor K=2:
γ₁(τ₁,τ₂,τ₃) = (2/3, 2/3, 5/3), orbit weights (16,4,1)/63 ⟹ S₂ = 2·(16·⅔ + 4·⅔ + 1·⅚)/63 = 2·(32+8+5)/189 = **10/21** ✓.

| K | engine S_K = 2Σ4⁻ᵐγ_{K−1}(τ_m) | frozen | verdict |
|---|---|---|---|
| 2 | 10/21 | 10/21 | ✅ |
| 3 | 31370/67963 | 31370/67963 | ✅ |
| 4 | 143195649659456490/308468774477179141 | = | ✅ |
| 5 | 2490699741144…/5350418720142… | = | ✅ |
| 6 | 3073594394508…/6593308371642… | = | ✅ |

The engine ⟹ frozen at every level through K=6. **Walk-back #31 not incurred** — the normalization
C_k(m) = 3[γ_{k−1}(τ_m) − γ_{k−2}(τ_m)] and the telescoping S_K = 2Σ4⁻ᵐγ_{K−1}(τ_m) are exact. The γ route lives.

## R9-B — γ-TABLES + X_n WELD (cross-thread, ALGEBRAIC forced): **WELD HOLDS**
γ_n(0) = X_n **exactly** for n = 1…5, byte-compared against the qx+1 corpus X_n = 1 + Σ_{j≤n} S_j:

| n | γ_n(0) = 3ⁿ‖μ_n‖² | corpus X_n = 1+Σ_{j≤n}S_j | weld |
|---|---|---|---|
| 1 | 5/3 | 5/3 | ✅ |
| 2 | 15/7 | 15/7 | ✅ |
| 3 | 177005/67963 | 177005/67963 | ✅ |
| 4 | 6626070796594781675/2159281421340253987 | = | ✅ |
| 5 | 18909241984277…/5350418720142… | = | ✅ |

**The τ=0 line of the collision profile IS the qx+1 corpus accumulation** — walk-back #27's A₁=5/3 sits on it, and
the whole X-column welds. **DC self-similarity confirmed:** the near-DC channel reads γ on the τ=0 line, so
C_k(DC) = 3[X_{k−1}−X_{k−2}] = **3·S_{k−1}** — the level-k ledger's DC channel is the level-(k−1) summit ×3
(C₂(DC)=3·⅔=2 ✓, C₃(DC)=3·10/21=10/7 ✓). The constant contains itself one level down, at Mersenne-collapsing
weight — "rationality from infinity" is now a literal identity.

**γ-table structure** (by stratum j=v₃(m), distinct values per stratum; DC = j=n is the τ=0 line):
- n=1: [j=0: 2/3] [DC: 5/3]
- n=2: [j=0: 24/49, 34/49, 40/49] [j=1: 10/7] [DC: 15/7]
- n=3: [j=0: 9 vals ~0.427..0.879] [j=1: 3 vals ~1.33..1.53] [j=2: 129950/67963] [DC: 177005/67963]
- n=4: [j=0: 27 vals ~0.405..0.917] [j=1: 9 vals] [j=2: 3 vals] [j=3: 1 val] [DC]
- n=5: [j=0: 81 vals ~0.401..0.953] [j=1: 27] [j=2: 9] [j=3: 3] [j=4: 1] [DC]
Within stratum j, γ takes exactly 3^{max(0,n−1−j)} distinct values (one per sub-orbit-class); the top bulk stratum
j=n−1 and the DC j=n are each a single value. Raw material for the pen's C̄_∞(j).

## R9-C — OFF-DC CONVERGENCE (measurement, NO fit, NO verdict)
γ_n(τ_m) for fixed m as n runs 1…5, verbatim (exact rationals computed; floats shown):

| m | v₃(m) | n=1 | n=2 | n=3 | n=4 | n=5 |
|---|---|---|---|---|---|---|
| 1 | 0 | 0.667 | 0.694 | 0.703 | 0.707 | 0.711 |
| 2 | 0 | 0.667 | 0.490 | 0.495 | 0.494 | 0.494 |
| 3 | 1 | 1.667 | 1.429 | 1.335 | 1.325 | 1.281 |
| 4 | 0 | 0.667 | 0.816 | 0.879 | 0.862 | 0.861 |
| 9 | 2 | 1.667 | 2.143 | 1.912 | 1.943 | 1.970 |

**Every fixed-m column is bounded** — none grows like X_n (which diverges 1.67→3.53 over the same range). The
higher-stratum columns (m=3, m=9) sit higher (~1.3, ~1.95) but stay bounded and settle. Consistent with Wilson's
pre-registered shape (convergence for each fixed m; divergence confined to τ=0); no column kills the re-scope. No
verdict drawn — pen adjudicates the limits.

## R9-D — NUMERAL WELDS (measurement, table only, NO conclusion)
**(i)** R66 primitive object P_n = Σ_{ξ prim mod 3ⁿ}|μ̂_n(ξ)|² (= Σ_{a,a′}μμ′c_{3ⁿ}(a−a′), computed directly)
vs this probe's ⟨γ_n⟩ = Σ_m 4⁻ᵐγ_n(τ_m)/(1/3):

| n | P_n = Σ_prim\|μ̂\|² | S_n (frozen) | ⟨γ_n⟩ | 3·S_{n+1}/2 | identities |
|---|---|---|---|---|---|
| 2 | 10/21 | 10/21 | 47055/67963 | 47055/67963 | P_n=S_n ✅; ⟨γ_n⟩=3S_{n+1}/2 ✅ |
| 3 | 31370/67963 | = | 214793474489184735/… | = | ✅ ✅ |
| 4 | 143195649659…/… | = | 3736049611716…/… | = | ✅ ✅ |
| 5 | 24906997411…/… | = | 4610391591763…/… | = | ✅ ✅ |

**Exact identifications (numerals, no verdict):** P_n = **S_n** level-by-level (R66's primitive object is the shell
itself); ⟨γ_n⟩ = **3·S_{n+1}/2** level-by-level (the next shell, scaled by the A∞ factor 3/2). They are distinct
objects at each level (P_2=10/21 vs ⟨γ_2⟩=47055/67963) and carry distinct limits (S_n→7/15; ⟨γ_n⟩→3·(7/15)/2=**7/10**,
= the A∞ mean twisted-collision density). Pen adjudicates whether "same object" or "shared numeral."

**(ii)** #27 chain amplitude: **NOT RECOVERABLE** from the archive — no frozen definition found (grep empty across
results/probes/STATE). Reported as NOT RECOVERABLE per spec; not reconstructed.

## Status
**R9 FULL PASS.** The collision identity S_K = 2Σ4⁻ᵐγ_{K−1}(τ_m) is verified exact from the μ tables, K=2…6
(R9-A gate); the τ=0 line welds to the qx+1 corpus X-column, n=1…5 (R9-B), confirming DC self-similarity
C_k(DC)=3S_{k−1}; every off-DC fixed-m column is bounded (R9-C, divergence confined to τ=0 as pre-registered);
and R9-D pins P_n=S_n, ⟨γ_n⟩=3S_{n+1}/2 exactly (→7/10), with #27 honestly NOT RECOVERABLE. The whole 7/15
campaign is now **one stationary function γ_∞ on ℤ₃**, evaluated on the dense ⟨4⟩-orbit {τ_m} with the
k-independent stratum weights W_j (R8-D3). **Still owed (pen):** γ_∞(τ) closed form and Σ_m 4⁻ᵐγ_∞(τ_m) = 7/30
(⟺ mean density 7/10 ⟺ S_∞ = 7/15). No fitting; exact rationals; nothing smoothed.
