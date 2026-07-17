# Result 9+10 (qx+1 paper) — PILLAR 3 REWRITTEN: `δ_q ≈ 0.82/ord_q(2)` is REFUTED (by 2.55×10¹³). The law is `δ_q = 2^{1−ord_q(2)}·(q−3)/q`.

**Date:** 2026-07-15. **Verdicts: R9 H_GATE CONFIRMED / R9 H_BASE CONFIRMED (base = 2.017 vs predicted 2) / R9 H_EXP INCONCLUSIVE (my rule was too weak) → R10 H_A_DEAD CONFIRMED / R10 H_CONCLUSIVE: B WINS OUTRIGHT.**

**Headline: pillar 3 goes from an 8-point empirical fit (R²=0.94, wrong by 10¹³) to a mechanism-derived closed form accurate to <1% at large q. And its `(q−3)/q` factor is pillar 2's constant — pillars 2 and 3 are not independent.**

Probes: `probe_9_pillar3_bakeoff.py`, `probe_10_pillar3_conclusive.py`. Logs: `result_9_bakeoff_log.txt`, `result_10_conclusive_log.txt`. Runtime: R10 = **1.6 s**.

## The mechanism (from R8, predicted BEFORE fitting)

R8 proved `π_k` is the q-adic self-similar measure and that **`ratio_1 = 0` exactly** — at k=1 the coding `v_1 ↦ 2^{−v_1} mod q` with `v_1 ∈ {1..ord_q(2)}` is a **bijection onto ⟨2⟩**. So:

> At k=2 a collision needs `2^{−v_2} ≡ 2^{−v'_2} (mod q)` ⟺ `v'_2 = v_2 + j·ord_q(2)`, `j ≥ 1`.
> Pair weight `~ 2^{−2v_2}·2^{−j·ord}`; summing `j ≥ 1`:
> **the cheapest collision costs a full period shift ⇒ overlap `~ 2^{−ord_q(2)}`.**
> `offdiag` sums over ORDERED pairs ⇒ every collision counted twice ⇒ prefactor **2**.

`1/ord_q(2)` has **no mechanism** — nothing in the structure produces a reciprocal.

## R9 — the bake-off, and the R² trap reproduced to three digits

**H_GATE CONFIRMED:** recomputed all 8 published δ from scratch (`c̃_q(2)=S_2/(q/3)²`, `S_2=X_2−X_1`) to ≤3.8e−4.

**H_BASE CONFIRMED — the sharp test.** Free 2-parameter fit `ratio_2 = a·c^{−ord}` recovered **`c = 2.01704`** against a mechanism that named **2** before fitting.

```
ratio_2 / 2^-ord :  2.76 2.66 2.12 2.77 2.20 2.08 2.00 2.02   <- spread 1.4x
ratio_2 * ord    :  0.69 1.00 0.02 0.01 0.07 0.32 0.04 0.11   <- spread 123x
```

**⚠️ H_EXP INCONCLUSIVE — second buggy decision rule, recorded not retro-fixed.** Log-space R²: exponential 0.99586 vs power 0.98022, Δ=0.0156 < my 0.05 threshold ⇒ INCONCLUSIVE by the pre-reg, and honoured as such. **The rule was badly designed:** ΔR² between two *free* 2-parameter fits on 8 monotone points across 165× cannot separate them. *(This is the second rule failure in two probes — H_DOM's "step<1.5" could not see linear growth. The fix, applied in R10: pre-commit NUMBERS, not adjectives.)* Note the competing power law has exponent **−4.4985, not −1** — it is not pillar 3's law; it is merely the thing that tied.

**★ H_R2TRAP — the methodological finding.** Recomputed linear fit:
```
LINEAR:  delta = -0.08903 + 0.81886*(1/ord)   linear R^2 = 0.94045
result_4 reported:            slope 0.82,             R^2 = 0.94
```
**Their slope and R² reproduced to 3 digits — while the same law is off by 121.5× at q=13, 53.4× at q=11, 24.3× at q=73.** Linear R² is dominated by the two largest points; the six small ones cluster near the origin where `1/ord` is also small. **R²=0.94 was real and meaningless.**

**★★ AND THE FITTED LAW CONTRADICTS A THEOREM.** Intercept `−0.089` ⇒ predicts **δ < 0 for ord > 9.2**. But R8's **H_LB proved δ > 0 always** (Cauchy–Schwarz per fiber). So `1/ord` dies both ways: with the intercept it violates positivity; forced through the origin it is off by 121×.

**Why the original "OOS validation" missed this:** `result_4`'s OOS primes q=31, 127, 73 have ord **5, 7, 9** — *interpolation inside* the fitted range (3–12). They tested the fit; they never tested the functional form.

## R10 — conclusive. Predictions pre-committed before the run.

Laws fit on the OLD primes (ord 3–12) only; new primes at ord **20, 23, 58** — far outside the fitted range. **No fitting on the new primes.**

**Self-checks (both PASS):**
- **Independent implementation:** 672,400-address combinatorial enumeration vs power iteration at q=41 → **2.776e−17**. Two unrelated routes to π.
- **Exact rational vs float:** `X_1` exact = float to **0.000e+00**.

| q | ord | measured `ratio_2` | **B** (`2.450·2^−ord`) | **P** (`73.87·ord^−4.4985`) | **A** (`0.82/ord`) |
|---|---|---|---|---|---|
| 41 | 20 | 1.9199e−06 | **1.22×** | 54.1× | ~2×10⁴× off |
| 47 | 23 | 2.4011e−07 | **1.22×** | 230.3× | ~10⁵× off |
| 59 | 58 | **−2.5e−16 = zero** | 8.5e−18 (sub-eps, as pre-stated) | predicted 8.6e−07 → **REFUTED** | predicted 0.01414 → **REFUTED** |

- **H_A_DEAD CONFIRMED — law A refuted by 2.55×10¹³.**
- **Law P REFUTED** — it predicted 8.6e−07 at q=59, nine orders above float noise; measured zero.
- **H_CONCLUSIVE: B WINS OUTRIGHT** — max |log10 miss| 0.085 (1.22×) vs P's 54–230×.

## ★★★ The 1.22× is a pure prefactor — and the prefactor is 2

B's miss is **1.22× at both** q=41 and q=47, identical to 3 digits. That is not a functional-form error. Backing the prefactor out of the new (large-q, uncontaminated) primes:

```
q=41: 1.9199e-6 / 2^-20 = 2.0132
q=47: 2.4011e-7 / 2^-23 = 2.0142
```

**The prefactor is 2** — exactly what the ordered-pair double-count predicts. The fitted `a = 2.450` was contaminated by the small-q primes (which scatter 2.1–2.8 from `O(1/q)` effects).

> ## `ratio_2 = 2^{1−ord_q(2)}` and **`δ_q = 2^{1−ord_q(2)}·(q−3)/q`**

`δ(2) = ratio_2·(q−3)/q` verified exactly at both new primes (0.927 = 38/41; 0.936 = 44/47). Checked against the OLD primes, which never informed it:

| q | ord | predicted | measured | miss |
|---|---|---|---|---|
| 73 | 9 | 0.003746 | 0.00375 | **0.1%** |
| 41 | 20 | 1.768e−6 | 1.780e−6 | **0.7%** |
| 47 | 23 | 2.232e−7 | 2.249e−7 | **0.8%** |
| 127 | 7 | 0.015256 | 0.01538 | **0.8%** |
| 31 | 5 | 0.05645 | 0.05871 | 4% |
| 11 | 10 | 1.420e−3 | 1.54e−3 | 8% |
| 17 | 8 | 0.006434 | 0.00722 | 12% |
| 7 | 3 | 0.1429 | 0.2103 | 47% |
| 13 | 12 | 3.756e−4 | 5.6e−4 | 49% |
| 5 | 4 | 0.05 | 0.0922 | 84% |

**Sub-1% for large q, degrading at small q exactly as an `O(1/q)` correction should.**

## What this means for the paper

- **Pillar 3 is rewritten**: `δ_q ≈ 0.82/ord_q(2)` (empirical fit, R²=0.94, wrong by 10¹³) → **`δ_q = 2^{1−ord_q(2)}·(q−3)/q`** (mechanism-derived base 2, matched prefactor 2, <1% at large q).
- **Pillars 2 and 3 are not independent** — pillar 3 carries pillar 2's `(q−3)/q` factor. With R7 (pillars 1↔2 via Pythagoras) and R8 (the "3" = `Σ_v p_v²`), **all three pillars now flow from the one self-similar-overlap mechanism.**
- **Still owed (flagged, not claimed):** the prefactor 2 is *matched* and has a plausible ordered-pair reason, but only the **cheapest** collision was derived — the full collision count is not done. The `O(1/q)` correction is unexplained. Both are counts, not fits, so both are gettable.

## Not at stake
THEOREM_C_745, Th 78.1–78.3, R81b, ε_k, R5's rate (pillar 1), R7's object identification, R6's blocked-route verdict.

_Reporting discipline: R9's H_EXP was INCONCLUSIVE by my own pre-registered rule and is reported as inconclusive, not laundered into a win — the rule's design flaw is named. R10 pre-committed exact NUMBERS for three unseen primes before running, so neither prediction could be moved afterward. The q=59 asymmetry (B's prediction is sub-machine-epsilon and cannot be confirmed there, only A and P refuted) was stated up front in the pre-reg, not discovered in the writeup. Two self-checks against independent implementations ran before the test. Author's structural priors this arc: now 5-for-10._
