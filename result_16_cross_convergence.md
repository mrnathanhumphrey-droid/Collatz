# Result 16 (qx+1 paper) — the cross-cell term CONVERGES at q≥5 and DIVERGES at q=3. That dichotomy IS domination. ⚠️ But my `3/q` rate prediction FAILED.

**Date:** 2026-07-15. **Verdicts: ★ H_CRIT ✓ CONFIRMED (the control: `D_{k+1}/D_k → 1.00092` at q=3) / ★ H_CONV ✓ CONFIRMED (cross converges at q≥5; Aitken limits reproduce R7's `C_q` to 0.4–0.9%) / ⚠️ H_RATE FAILED — my prior lost.**

**Headline: the domination DICHOTOMY is now measured directly — cross converges for q≥5, diverges linearly at q=3. But the *rate* is not `3/q`, so there is still no quantitative bound.**

Probe: `probe_16_cross_convergence.py`. Log: `result_16_cross_convergence_log.txt`. Runtime: **12.6 s**, peak **1.66 GB** (dense chain — the transition matrix is ~100% dense, so dense is *cheaper* than sparse here).

## ★ H_CRIT — the control, and it is rock solid

At q=3, successive-difference ratios `D_{k+1}/D_k` (where `D_k := ratio_{k+1} − ratio_k`; if `deficit ~ A·r^k` then `D_{k+1}/D_k → r`):

```
0.78510 → 1.00267 → 1.00140 → 0.99855 → 1.00092
```

**Dead on 1.0 from k=4 — off by 0.1%.** No decay, ever ⇒ linear divergence. The *same statistic* that must read `<1` at q≥5 reads exactly `1` at q=3. The control passes, so the framing is sound.

## ★ H_CONV — confirmed, and it independently validates R7

Cross **converges** at q≥5 (sequences increasing, finite Aitken limits). Against R7's `C_q − 1 = c̃_q·q/(q−3) − 1`, computed by an entirely different route:

| q | Aitken limit | `C_q − 1` (R7) | off | cross limit |
|---|---|---|---|---|
| 5 | 0.221273 | 0.2193 | 0.9% | 0.087938 |
| 7 | 0.366103 | 0.3695 | 0.9% | 0.080387 |
| 11 | 0.002112 | 0.00212 | **0.4%** | 0.000157 |

> **The DICHOTOMY is confirmed: cross converges for q≥5, diverges linearly at q=3. That dichotomy IS domination** — and R7's identity gets a free independent check.

## ⚠️ H_RATE — FAILED. My prior lost.

| | q=5 | q=7 | q=11 |
|---|---|---|---|
| observed `D_{k+1}/D_k` | 0.630 | 0.392 | 0.527 |
| predicted `3/q` | 0.600 | 0.429 | 0.273 |
| verdict (pre-registered rule) | INCONCLUSIVE (5% off but drifting **away**) | CORROBORATED (8.6%, toward) | INCONCLUSIVE (**93% off**) |

**The observed rates are not even monotone in q** — 0.63, 0.39, 0.53 — which no `3/q` law can produce.

**Why I was wrong (diagnosed, not excused):** my derivation assumed `M_j/(q/3)^j = c̃_q` **exactly** for every j. It isn't — `c̃_q(j)` has its own convergence toward `c̃_q`, and that contaminates the deficit. The observed decay is `max(3/q, rate of c̃_q(j))`, not `3/q`. I derived the geometric tail and forgot the coefficient was still moving.

## Where Result 1 stands

**Confirmed:** within-cell bounded + closed-form (R15); cross-cell converges at q≥5, diverges at exactly q=3 (this). **That dichotomy is what domination needs, and it is now measured rather than inferred.**

**NOT established:** the *rate* of convergence ⇒ **no quantitative bound** ⇒ **Result 1 is still not a theorem.** A proof needs a genuine bound on the cross-cell **character sum** (R13's object), not a rate guess.

**Sixth independent sighting of the phase boundary**, and the sharpest yet: `D_{k+1}/D_k → 1` at q=3 versus `< 1` at q≥5.

## Pattern worth recording

My **structural** priors keep landing (the grading; the tower closed form; the cross/within split; the q=3 slope 7/15 to 0.06%). My **quantitative** priors keep losing (the replacement prefactor's reason; `O(1/q)`; "cross stays flat"; `3/q`). Five mis-specified decision rules. Consistent enough to be a rule: *see the structure, don't guess the numbers on it* — which is exactly why pre-committing numbers keeps paying.

## Not at stake
R10's law, R11, R13, R14, R15's within-cell identity, R5's rate, R6, R7, R12, THEOREM_C_745, Th 78.1–78.3, R81b, ε_k.

_Reporting discipline: rates were committed as NUMBERS before the run. The threshold was deliberately loose (20%) and labelled as a SIGN CHECK on the mechanism rather than a measurement — with 2–3 difference-ratios a tight threshold would be theatre, and that was said in the pre-registration, not afterward. H_RATE fired INCONCLUSIVE/CORROBORATED/INCONCLUSIVE and is reported as a FAILURE of my prediction rather than salvaged from the one q that corroborated. H_CRIT was included specifically as a control so that a broken statistic would be caught. Resource use stated up front. Author's structural priors this arc: 17-for-26._
