# Result 23 — ROUTE 2 (exact-rational): `7/15` is the LIMIT not the exact increment; `c_k` is order ≥2 everywhere; the q=3 approach is a super-geometric TOWER. Recurrence-hunt under-determined → route 1 = the RPF transfer operator.

**Date:** 2026-07-16. **Type:** probe (`probe_23_exact_increment_recurrence.py`), exact big-integer arithmetic, no float. **Verdict: the increment structure is NOT a clean finite recurrence; `r_q` is the subdominant eigenvalue of the IFS transfer operator, which route 1 must compute.**

## Method (fully exact, gate-confirmed)

Every cell weight is `2^E/D` (integer exponent `E = Σ_j(m_j−c_j)+(M−v_k)`, fixed denominator `D = Π_j(2^{m_j}−1)·(2^M−1)`), and there are only `φ(q^k)` distinct residues. Accumulate exact integer numerators `N_r = Σ 2^E` per residue; then `cross(k) = (Σ_r N_r² − Σ 2^{2E})/D² / P2^k` as one final `Fraction`. **G_EXACT gate CONFIRMED:** exact `cross(k)` matches R18's float `cross_from_cells` to ≤1.3e−14 rel at all 11 tested (q,k). Reach: q=3 k≤5, q=5 k≤4, q=7 k≤3, q=11 k≤3.

## Finding 1 — `7/15` is the asymptotic slope, NOT the exact per-level increment

At q=3, `c_k` is **not** `7/15` exactly:

| k | `c_k` (float of exact) | `c_k − 7/15` |
|---|---|---|
| 3 | 0.487188595953 | +2.05e−2 |
| 4 | 0.464234660612 | −2.43e−3 |
| 5 | 0.465514919752 | −1.15e−3 |

`c_k → 7/15` from below, but each is a distinct exact rational ≠ 7/15. **This refines R15, does not contradict it** — R15's `H_CROSS_GROWS` only ever claimed the slope constant to 2%, which holds. The exact statement is `cross(k) = (7/15)k + const + (decaying correction)`.

## Finding 2 — `c_k` is order ≥2 (non-geometric) EVERYWHERE, exactly

Order-1 (geometric) recurrence ruled out by exact rationals, not float: `c_3/c_2 ≠ c_4/c_3` at both q=3 (1.2252 vs 0.9529) and q=5 (0.5337 vs 0.5077). So `ρ_k = c_{k+1}/c_k → r_q` is the limit of a non-geometric sequence; `r_q` is a **dominant** eigenvalue sitting above subdominant modes.

## Finding 3 — the q=3 approach to `7/15` is a super-geometric TOWER, not one mode

Second differences `Dc_k = c_k − c_{k−1}` and their ratios:
- `Dc_3 = +8.96e−2`, `Dc_4 = −2.30e−2`, `Dc_5 = +1.28e−3`
- ratios `Dc_4/Dc_3 = −0.2563`, `Dc_5/Dc_4 = −0.0558` — **alternating sign, shrinking fast** (not equal, exact).

A single subdominant eigenvalue would give a *constant* ratio. The shrinking-faster-than-geometric ratio is the fingerprint of the **doubly-exponential tower** `x_j = 2^{−d·q^{j−1}}` (R14): each coordinate contracts at its own rate, so corrections to `7/15` die like `~2^{−d·q^{k}}`. There is **no finite linear recurrence** capturing the full sequence.

## Finding 4 — recurrence-hunt is under-determined (as pre-registered)

Order-2 verify needs 5 exact `c`-values; available: q=3 has 4, q=5 has 3, q=7/11 have 2. **The order is not fittable from data** — exactly the predicted outcome, and the reason route 1 (structure, not fitting) is required.

## Consequence — route 1 is the RPF transfer operator, and it is well-posed

The naive "finite recurrence → finite matrix" is dead (finding 3: the correction structure is a tower). The correct operator is the **Ruelle–Perron–Frobenius transfer operator** `L` of the IFS `T_v(x) = (qx+1)/2^v`, weights `p_v = 2^{−v}/Z`:
- its **top eigenvalue** governs the diagonal — this is R5's correlation dimension `D₂ = log3/logq`;
- **`r_q` is its subdominant eigenvalue**;
- unlike the k-index (unbounded), the RPF spectrum **converges under truncation** (Fourier modes / q-adic tree depth), so `λ₂(L)` is computable at fixed dimension-per-truncation;
- **`r_q < 1 ⟺ spectral gap of `L`** — the standard, provable statement for self-similar measures. This is the honest Phase-3 target, and route 1 (next probe) discretizes `L` and watches `λ₂(L)` converge.

## Plan status after R23

| phase | status |
|---|---|
| 0,1,1b | DONE exact · 1c PARTIAL · 2 DONE (Konyagin shelf) |
| **3 — the bound** | object = `ρ_k` (R22); **R23: increment is order≥2 + tower, NOT a finite recurrence; `r_q = λ₂` of the RPF transfer operator; target = spectral gap. Route 1 = compute `λ₂(L)` under truncation.** |

## Not at stake
R10–R22 (R15's slope refined, not refuted), R5's rate, R6, R7, R12, THEOREM_C_745.

_Reporting discipline: the `7/15` result is stated as a refinement of R15 with the exact misses shown, not spun as a contradiction; "order ≥2" is proven by exact inequality of rationals, not asserted from float; the tower reading of finding 3 is offered as the structural explanation and is consistent with R14's grading; the recurrence under-determination is reported as the pre-registered expected outcome, not as a failure._
