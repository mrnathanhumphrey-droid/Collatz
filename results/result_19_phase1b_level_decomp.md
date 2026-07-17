# Result 19 — PHASE 1b: ★ the level structure. `cross(k) = cross(k−1) + c_k`, and `c_k` is the `j≥1` mass — which at q=3 IS R7's `M_k`.

**Date:** 2026-07-15. **Verdicts: H_GATE ✓ (classifier exact) / ★ H_NEWLEVEL — structure FOUND (reported as measurement; no verdict, no fit, per pre-reg) / ⚠️ q=3,k=3 does NOT fit.**

**Headline: the target "bound cross(k) in k" becomes "bound `Σ_k c_k`", one subgroup count per level — and `c_k` is exactly the mass of pairs whose finest coordinate differs (⟺ R13's `j ≥ 1`).**

Probe: `probe_19_phase1b_level_decomp.py`. Log: `result_19_phase1b_log.txt`. Runtime: **1.4 s**.

## ★ H_NEWLEVEL — the level structure

| q | k | cross(k) | no-new-coord | has-new-coord | cross(k−1) | match |
|---|---|---|---|---|---|---|
| 3 | 4 | 1.34905639 | **0.88484198** | 0.46421441 | **0.88482173** | 2e−5 |
| 5 | 3 | 0.05962378 | **0.03887536** | 0.02074841 | **0.03887521** | 4e−6 |
| 7 | 3 | 0.06804208 | **0.04702543** | 0.02101665 | **0.04702534** | 2e−6 |
| **3** | **3** | 0.88482173 | **0.42327396** | 0.46154777 | **0.39763314** | **6% — DOES NOT FIT** |

Pairs that do **not** touch the new finest coordinate reproduce `cross(k−1)`; everything new is carried by pairs where `v_k` differs. And "differs in `v_k`" ⟺ **`j ≥ 1`** (the shift is a multiple of `d`) — so **R13's shift index IS the level index.** Cross-check: the `j=0` share at q=5,k=3 is 0.652, and `0.348 × cross(3) = 0.0207` = the has-new-coord mass exactly.

> **`c_k` := the `j ≥ 1` mass at level k. &nbsp; `cross(k) = cross(k−1) + c_k`.**

## `c_k` is R7's object

- **q=3:** `c_3 = 0.4615`, `c_4 = 0.4642` — essentially constant at ≈ **7/15** ⇒ **linear divergence**, matching R15's slope 0.46577.
- **q=5:** `c_3 = 0.0207` — *precisely* R16's measured increment sequence (0.0207, 0.0106, 0.0066, 0.0042).

At q=3, `total(k) = X_k − 1`, so `c_k = X_k − X_{k−1} = M_k(1)` — **R7's identity.** The level increment of the cross term *is* the paper's central object, arriving from a fourth independent direction.

## ⚠️ The anomaly, reported not buried

**q=3, k=3 misses by 6%** (0.4233 vs 0.3976) while q=3, k=4 matches to 2e−5. Possibly a genuine transient — q=3 has `d=2`, and R15 independently found k=3 anomalous there (`x_3 = 2^{−18}` is not negligible at d=2). Or the split is not clean at small k. **No verdict claimed:** the pre-registration committed to reporting H_NEWLEVEL as a measurement with no fit, precisely so an anomaly like this cannot be smoothed over.

## Supporting structure (measurements, no verdicts)

- Mass concentrates on pairs differing in **exactly 2 coordinates**: 0.74–0.90 share.
- **`v_1`** (the coarsest) is the lowest differing coordinate for 0.91 of the mass at q=5,k=3.
- Consistent with the geometric weights — but "consistent with" is not evidence.

## Where the plan stands

| phase | status |
|---|---|
| 0 — general-k cross expression | **DONE**, exact |
| 1 — k=2 subgroup form | **DONE**, exact |
| **1b — level structure** | **DONE**: `cross(k) = cross(k−1) + c_k`, `c_k` = the `j≥1` mass |
| **1c — REQUIRED next** | R13's condition is **k=2-specific**; `c_k` at general k needs the analogous condition from the k-term value. **A derivation, not a probe.** |
| 2 — literature | pending; target confirmed correct |
| 3 — the bound `c_k ≤ C_q·r_q^k` | the math |

## Not at stake
R10, R11, R13, R14, R15, R16, R17, R18. Exploratory structure-finding.

_Reporting discipline: only the GATE carried a threshold, and it was an exact identity — deliberate, since six of my pre-registered thresholds have been mis-specified this arc, every one blind to a parameter of the mechanism. H_NEWLEVEL/H_SHALLOW/H_Q3 were pre-committed as measurements with NO verdict and NO fit; the q=3/k=3 miss is reported at full weight rather than averaged away. Budget skips were announced, not silent._
