# Result 22 — PHASE 3 (reconciliation): the weighted increment ratio `ρ_k = c_{k+1}/c_k` is the correct Phase-3 object. It → 1 at q=3 (divergence) and < 1 for q≥5 (convergence). The R20/R21 caveat is RESOLVED. My committed value `3/q` is REFUTED.

**Date:** 2026-07-16. **Type:** probe (`probe_22_phase3_increment_ratio.py`). **Verdict: the object is pinned and the 1/q-vs-divergence paradox is gone; the *value* `r_q` for q≥5 is not (and need not be) determined by computation.**

## Headline

The quantity that governs `cross(k)` is **not** R20's unweighted per-level collision rate (~1/q). It is the **weighted increment ratio** `ρ_k := c_{k+1}/c_k`, where `c_k = cross(k) − cross(k−1)` (R19's new-finest-coordinate mass). Measured via the stationary route (R18 H_P0 gate: `cross(k) = ‖π_k‖²/P2^k − 1 − ratio_within(q,k)`):

| q | `ρ_k` sequence | behavior | `cross` |
|---|---|---|---|
| **3** | 1.225, 0.953, **1.003, 1.001, 0.999, 1.001** | `→ 1.000` (c_k locks on 7/15) | linear divergence ✓ |
| 5 | 0.534, 0.508, 0.624 | rising, < 1 | converges ✓ |
| 7 | 0.447, 0.357 | < 1 | converges ✓ |
| 11 | 0.181 | < 1 | converges ✓ |

## The caveat (R20/R21) is resolved

The flagged paradox: R20's per-level rate ~1/q would give `c_k ~ (1/q)^{k−1}`, which **converges at q=3 too** (1/3 < 1) — contradicting the q=3 divergence. Resolution: **1/q and `ρ_k` are different objects.** 1/q is the *unweighted* per-level collision probability; `ρ_k` is the *weighted* increment ratio. `ρ_k` is the one that controls `cross(k)`, and at q=3 it is **1.000, not 1/3** — which is exactly why q=3 diverges. `c_k` at q=3 locks onto **7/15** (R15's slope), giving `ρ_k → 1` cleanly (every k≥4 ratio within 0.3%, numerical wobble on the power-iteration norm). No contradiction remains. **Do-not-assume flag from R21 is retired: `c_k` is NOT `(1/q)^{k−1}`; it is `c_2·∏ρ_j` with `ρ_j → 1` at q=3.**

## G2 (phase boundary) — CONFIRMED

`ρ_k → 1` at q=3, `ρ_k < 1` at q=5, 7, 11, at every reachable k. The phase boundary shows up cleanly in this object. **The revised Phase-3 target is validated as the correct statement:** prove `ρ_k ≤ r_q < 1` uniformly in k for q≥5.

## PRED `ρ_q → 3/q` — REFUTED (8th consecutive quantitative-prior miss)

Committed before the run: `ρ_q = 3/q` (would reconcile: 1 at q=3, <1 for q≥5, repurposing R16's dead 3/q). It missed with **inconsistent signs**: q=5 overshoots (0.624, still *rising* past 0.60), q=7 undershoots (0.357 vs 0.429), q=11 undershoots (0.181 vs 0.273). No clean law. Structural prior (G2) held; quantitative prior lost — the invariant pattern of this arc (structural ~20/many, quantitative ~0/8). **The reconciliation never depended on the value**, so this miss costs nothing.

## The honest open piece + the one flag

- **`r_q` for q≥5 is not pinned** at reachable k (local budget n·M ≤ 25e6 → q=5 caps at k=5). The sequences are pre-asymptotic.
- **⚠️ FLAG: the q=5 ratio is CLIMBING** (0.51 → 0.62 → 0.64). `ρ<1` at reachable k is *consistent* with convergence but does **not** secure the uniform bound `r_q<1` empirically — a rising sequence could in principle approach 1. This is not a probe defect; it is precisely the open analytic content of Phase 3. The uniform-in-k bound is a **spectral-gap statement about a fixed operator**, provable at fixed dimension — not something finite-k computation can close.

## Numerical hygiene

`c_k = cross(k) − cross(k−1)` is a difference of floats; guarded to flag `|c_k| < 1e-9` as cancellation noise. In all k-ranges reported, `c_k > 1e-3` (q=3..11 down to 2e-5 at q=11 k=3, still 4 sig-figs clean). `cross(1) = 0` exactly at every q (no cross-cell at k=1), so `c_2 = cross(2)`.

## Plan status after R22

| phase | status |
|---|---|
| 0, 1, 1b | DONE, exact |
| 1c | PARTIAL (level 2 verified at scale; levels ≥3 not derived) |
| 2 — literature | DONE (Konyagin shelf, wrong axis) |
| **3 — the bound** | **object PINNED (`ρ_k`), caveat RESOLVED, G2 confirmed; the uniform `r_q<1` bound is now a clean spectral-gap target, value TBD by proof not compute** |

## Not at stake
R10–R21, R5's rate, R6, R7, R12, THEOREM_C_745. This measures the Phase-3 object and retires a caveat; it unbanks nothing.

_Reporting discipline: the phase boundary (G2) is the confirmed structural result and is stated as such; the `3/q` miss is reported as a miss with its inconsistent signs, not spun; the rising q=5 sequence is flagged as an open risk to the uniform bound rather than buried, because "ρ<1 at k≤5" is exactly the kind of finite-k pass that must not be oversold as the theorem._
