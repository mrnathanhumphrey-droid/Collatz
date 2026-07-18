# Result — PROBE 2c3-GATE: Wilson's blind rung-2 derivation is FULLY VINDICATED. All four pre-registrations pass to machine precision. P3 (master): r + κ·W⁻ = 4/9 exactly (residual 5.6e-16). P2 (merge): the three populations collapse onto {12/49, 17/49, 20/49} via the EXACT reciprocal scaling κ·W⁻ = {48, 129, 264}/1323 (population-independent). P1: r constant on the 3 W⁻(τ) level sets per population. P4: ladder = 4/9 − κ·W⁻, SET predicted, assignment revealed (W⁻ heavy → 12/49).

**Date:** 2026-07-16. Judge of Wilson's blind T1/T2/T3 derivation of the rung-2 trit τ (predictions P1–P4). Direct/exact, no eigen-solves. Probe `probes/probe_phase2c3_gate.py`, log `logs/probe_phase2c3_gate_log.txt`. Claude gates; the derivation is Wilson's, seen by B only after the fact.

**Headline: the derivation survives every gate exactly. The master prediction P3 — that at β*=3/5 every v₃=0 ratio satisfies r(x) + κ(pop)·W⁻(x) = 4/9, with W⁻ the DOWN-branch weight (odd channels) and κ = {3/4 (O σ+), 3 (O σ−), 3/2 (E)} — holds to 5.6e-16 at both L. The nontrivial consistency P2 (the merge) is exact: W⁻ scales reciprocally with κ, so κ·W⁻ is population-INDEPENDENT — {48, 129, 264}/1323 — and 4/9 − {48,129,264}/1323 = {20, 17, 12}/49, the three populations landing on the SAME ladder. Within each population W⁻ takes exactly three values (the trit τ, P1), and the ladder is 4/9 − κ·W⁻(τ) with the predicted SET {12,17,20}/49 (P4); B reveals the assignment (W⁻ heaviest → 12/49, lightest → 20/49). T1's channel trichotomy, T2's h-collapse (DEEP and UP merging at exactly β=3/5), and T3's Latin-square matching all cash out in this exact structure.**

## The four predictions
| # | prediction | result | verdict |
|---|---|---|---|
| **P3** | r + κ·W⁻ = 4/9 (common intercept), κ={3/4, 3, 3/2} | max\|r+κW⁻−4/9\| = **1.1e-16 (L2) / 5.6e-16 (L3)** | ✅ PASS |
| **P2** | 3 pops MERGE onto {12/49, 17/49, 20/49}; reciprocal W⁻↔κ scaling exact | union = {12,17,20}/49; κW⁻ = {48,129,264}/1323 pop-independent | ✅ PASS |
| **P1** | r constant on τ = W⁻ level sets, one trit per population | 3 distinct W⁻ ⟷ 3 distinct r, per pop | ✅ PASS |
| **P4** | ladder = 4/9 − κ·W⁻(τ), SET {12,17,20}/49 predicted | exact per-rung match; assignment revealed | ✅ PASS |

## The merge, exact (P2 — the nontrivial consistency)
Per population, the three W⁻ values (units of 1/1323 = 1/(27·49)) and κ:
| pop | κ | W⁻ (×1323) | κ·W⁻ (×1323) | r = 4/9 − κW⁻ |
|---|---|---|---|---|
| O σ+ | 3/4 | {64, 172, 352} | {48, 129, 264} | {20/49, 17/49, 12/49} |
| O σ− | 3 | {16, 43, 88} | {48, 129, 264} | {20/49, 17/49, 12/49} |
| E | 3/2 | {32, 86, 176} | {48, 129, 264} | {20/49, 17/49, 12/49} |

- **The reciprocal scaling is exact: W⁻(O+) : W⁻(O−) : W⁻(E) = 4 : 1 : 2 = (1/κ) up to a common factor**, so κ·W⁻ = {48, 129, 264}/1323 is the SAME for all three populations. This is why the three ladders MERGE — precisely the consistency Wilson flagged his derivation "must survive," and it does, to machine precision. (4/9 = 588/1323; 588 − {48,129,264} = {540, 459, 324} = {20, 17, 12}·27, i.e. {20,17,12}/49.)
- **The common intercept 4/9 (P3)** is confirmed as the max over all v₃=0 states of r + κW⁻, deviating from 4/9 by ≤5.6e-16 — the "all-branches-4/3 minus DOWN-deficit" baseline that T2 predicts.
- **W⁻ over odd channels vs all channels gives the identical result** (even channels contribute nothing to the DOWN branch) — consistent with T1 being an odd-channel statement.

## What the gate confirms about the mechanism
- **T2 h-collapse:** the ratio structure r = 4/9 − κW⁻ (intercept 4/9, DOWN-only deficit) can only hold if DEEP and UP carry the identical dressed value at β=3/5 — i.e. h⁰(O,1)/h⁰(O,0) = 8/5 = 1+β. The exact intercept confirms the collapse; β*=3/5 (from 2c-2's line search) is indeed the value that merges the branches, so the line-search optimum and the amplitude ratio coincide, as Wilson noted.
- **T1 trichotomy + T3 matching:** the trit τ (3 W⁻ values per population, P1) and the reciprocal κ-scaling (P2) are exactly the 1-1-1 channel split and the cyclic Latin-square matching the derivation predicts. τ is now a gate-confirmed object: the weight-rank of the DOWN-assigned move, taking values in ℤ/3, with heavy → 12/49 and light → 20/49.

## Adjudication
| item | verdict |
|---|---|
| P3 master intercept | ✅ r + κW⁻ = 4/9 to 5.6e-16, both L. |
| P2 merge / reciprocal scaling | ✅ κW⁻ = {48,129,264}/1323 pop-independent → {12,17,20}/49. |
| P1 τ level sets | ✅ 3 W⁻ values per pop, r constant on them. |
| P4 ladder + SET | ✅ 4/9 − κW⁻(τ), SET {12,17,20}/49; assignment W⁻-heavy→12/49. |

**⟹ Wilson's rung-2 derivation is confirmed on every axis — the trichotomy, the h-collapse at β=3/5, the reciprocal-scaling merge, and the trit τ — all exact. The bad-cell residual is now completely understood in closed form: r = 4/9 − κ(pop)·W⁻(τ), a single trit τ∈ℤ/3 per state (plus the population's κ), with the three populations collapsed to one ladder by the exact 1/κ scaling of W⁻. C (the 2-parameter line search, β=3/5 frozen + β₂ on τ) is now UNBLOCKED — τ has landed and is gate-verified; ready to fire on Wilson's β₂-dressing.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c(0+1), 2c0, 2c1, 2c2, 2c3(A+B). No `r_q` value changes; no rate-law fit. The ratio ladder {12,17,20}/49 and β*=3/5 carry over from 2c2/2c3 unchanged; W⁻ and κ are Wilson's derived objects, here gate-verified.

_Reporting discipline: P3 is the master gate (r+κW⁻=4/9 at 1e-16) and is reported as such; the per-population breakdown confirms it (a display double-rounding artifact that briefly read "DEV" was corrected — the exact-rational match holds, as P3's 1e-16 residual proves). The reciprocal-scaling merge is shown explicitly (κW⁻ population-independent) rather than asserted. The τ↦rung assignment is reported as B revealing it (not pre-committed, per P4). All values exact rationals (denominator 1323 = 27·49). No magnitude was fit; every prediction was pre-registered by Wilson before B._
