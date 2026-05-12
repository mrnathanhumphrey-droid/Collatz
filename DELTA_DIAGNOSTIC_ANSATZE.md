# DELTA_DIAGNOSTIC_ANSATZE — five pre-registered ansatz fits with held-out residual

**Date:** 2026-05-12. Wilson — leading-vs-subleading diagnostic, Phase 2 of 3.

## Protocol

- **Fit on n=2..5**; **n=6 held out** as the falsifier point.
- Five ansatze pre-specified in the brief; no others tried.
- "Fits" defined as |predicted_δ_6 − actual_δ_6| / |actual_δ_6| < 0.20 (20% threshold, pre-registered, looser than usual because n=2..5 is a 4-point sample).
- Inputs from QUANTITIES.md:
  - δ_2 = +0.00476190
  - δ_3 = +0.00740250
  - δ_4 = +0.00590310
  - δ_5 = +0.00352297
  - δ_6 = −0.00146744 (held out)
- The 4 fit points are all positive; the held-out point is **negative** — the sign-flip is the falsifier.

---

## Ansatz (a) — geometric: δ_n = c · ρ^n

**Method.** Least-squares fit of log|δ_n| = log|c| + n·log|ρ| over n=2..5 (sign treated as positive across fit range).

**Fit.** Linear regression on (n, log|δ_n|):
- slope ≈ −0.099 → ρ ≈ 0.906
- intercept ≈ −4.912 → c ≈ 0.00735

**Prediction at n=6.** δ_6^pred = 0.00735 · 0.906⁶ = 0.00735 · 0.5552 = **+0.00408**.

**Actual.** δ_6 = −0.00147.

**Residual.** |0.00408 − (−0.00147)|/0.00147 = 0.00555/0.00147 ≈ **3.78 (378%)**.

**Result:** WRONG SIGN, off by ~4×. **FAILS** the 20% threshold by an order of magnitude. **Inherently** cannot match a sign-flip at the held-out point with real ρ.

---

## Ansatz (b) — power-law correction: δ_n = c · n^(−α) · (1/2)^n

**Method.** Equivalent to δ_n·2^n = c·n^(−α). Take log: log(δ_n·2^n) = log c − α·log n. Linear regression on (log n, log(δ_n·2^n)).

**Fit values δ_n·2^n on n=2..5:** 0.01905, 0.05921, 0.09446, 0.11274 (positive, monotone increasing).

Linear regression in log-log:
- slope ≈ +3.87 → α ≈ **−3.87** (negative, i.e., the "correction" grows like n^{+3.87})
- log c ≈ −7.46 → c ≈ 5.74×10⁻⁴

**Prediction at n=6.** δ_6^pred = 5.74×10⁻⁴ · 6^{3.87} · (1/64) ≈ 5.74×10⁻⁴ · 1028 / 64 ≈ **+0.00922**.

**Residual.** |0.00922 − (−0.00147)|/0.00147 ≈ **7.27 (727%)**.

**Result:** WRONG SIGN, off by ~7×. **FAILS** badly. Note that α is negative — the data over n=2..5 wants δ_n·2^n to *grow*, not decay, in n. There is no exponent α that produces a sign flip at n=6.

---

## Ansatz (c) — log correction: δ_n = c · log(n) · (1/2)^n

**Method.** Linear in c. Closed-form LSQ: c = Σ(δ_n · g_n) / Σ(g_n²), where g_n = log(n)·(1/2)^n.

**Computation:**
- g_n for n=2..5: 0.1733, 0.1373, 0.0866, 0.0503
- Σ(δ_n·g_n) = 0.000825 + 0.001016 + 0.000511 + 0.000177 = 0.002529
- Σ(g_n²) = 0.05891
- **c ≈ 0.04293**

**Prediction at n=6.** δ_6^pred = 0.04293 · log(6) · (1/64) = 0.04293 · 1.7918 · 0.015625 ≈ **+0.001202**.

**Residual.** |0.001202 − (−0.00147)|/0.00147 ≈ **1.82 (182%)**.

**Result:** WRONG SIGN, off by ~2×. **FAILS.** Log-correction is strictly positive for n≥2 with this c sign, cannot match sign flip.

---

## Ansatz (d) — two-term superposition: δ_n = c_1·ρ_1^n + c_2·ρ_2^n

**Method.** Four parameters, four fit points (n=2..5) → exact fit. Use **Prony's method**: enforce δ_{n+2} = a·δ_{n+1} + b·δ_n linear recurrence; solve for (a,b) from n∈{2,3}; ρ_1, ρ_2 are roots of x² − a·x − b = 0; then back out (c_1, c_2) from initial conditions.

**Linear system:**
- δ_4 = a·δ_3 + b·δ_2  →  0.00590 = 0.00740·a + 0.00476·b
- δ_5 = a·δ_4 + b·δ_3  →  0.00352 = 0.00590·a + 0.00740·b

**Solve.** Determinant = 0.00740² − 0.00476·0.00590 = 0.00002668. 
- a ≈ +1.0082
- b ≈ −0.3284

**Characteristic polynomial:** x² − 1.0082·x + 0.3284 = 0.
**Discriminant:** 1.0082² − 4·0.3284 ≈ −0.2971 (**negative**).
**Roots:** complex conjugate pair ρ_{1,2} = 0.5041 ± 0.2725·i.
- |ρ| = √(0.5041² + 0.2725²) = √0.3284 ≈ **0.5731**
- arg(ρ) = arctan(0.2725/0.5041) ≈ **0.495 rad** (28.4°)

**Held-out prediction via recurrence:**
δ_6^pred = a·δ_5 + b·δ_4 = 1.0082·0.00352 − 0.3284·0.00590 = 0.003549 − 0.001938 = **+0.001611**.

**Residual.** |0.001611 − (−0.00147)|/0.00147 = 0.003081/0.00147 ≈ **2.10 (210%)**.

**Result:** WRONG SIGN. **FAILS.** Notable: Prony gives **complex-conjugate roots**, meaning the "two-term" fit is *intrinsically oscillating* (an instance of ansatz e in disguise, but with damping rate |ρ|=0.5731 instead of forced 1/2). The cos-period is 2π/0.495 ≈ 12.7 steps — only half a cosine cycle fits in n=2..6, so the held-out prediction lies just before the first zero crossing of the implied cosine, near its peak, hence positive.

The Prony fit reproduces n=4,5 exactly by construction. The failure to predict n=6 is the structural failure of the two-term form, not a numerical fit issue.

---

## Ansatz (e) — oscillating: δ_n = c · cos(ω·n + φ) · (1/2)^n

**Method.** Three parameters, four fit points. Use the cosine sum identity exactly on n=2,3,4:
cos(α_{n+1}) + cos(α_{n−1}) = 2·cos(ω)·cos(α_n), where α_n = ω·n + φ.

**Setup.** Let A_n = δ_n·2^n = c·cos(α_n). Values: A_2 = 0.01905, A_3 = 0.05921, A_4 = 0.09446, A_5 = 0.11274.

Using A_2 + A_4 = 2·cos(ω)·A_3:
- 0.01905 + 0.09446 = 2·cos(ω)·0.05921
- 0.11351 = 0.11842·cos(ω)
- **cos(ω) ≈ 0.9585** → ω ≈ ±0.288 rad

**Prediction at n=5 (consistency check):**
A_5^pred = 2·0.9585·A_4 − A_3 = 0.18105 − 0.05921 = 0.12184
δ_5^pred = 0.12184/32 = 0.003807 (actual 0.003523 — close, ~8% off)

**Prediction at n=6:**
A_6^pred = 2·0.9585·A_5 − A_4 = 2·0.9585·0.11274 − 0.09446 = 0.21612 − 0.09446 = 0.12166
δ_6^pred = 0.12166/64 ≈ **+0.001901**.

Alternative LSQ over n=2..5 (3 free parameters, slight overdetermination) gives a similar number — local minimum around (c ≈ 0.13, ω ≈ 0.29, φ ≈ −1.4) producing δ_6^pred ≈ +0.0019–0.0022.

**Residual.** |0.00190 − (−0.00147)|/0.00147 ≈ **2.29 (229%)**.

**Result:** WRONG SIGN. **FAILS.** The cosine recurrence with cos(ω) ≈ 0.96 forces δ_6 to remain in the *first half-period* (still near its peak), not yet at zero crossing. To force a zero crossing at n=5.5, ω would need to be roughly π/4 ≈ 0.785 rad, but that's inconsistent with the n=2..5 progression (which constrains cos(ω) ≈ 0.96, i.e., ω small).

---

## Summary table

| Ansatz | Fit params | Pred δ_6 | Actual δ_6 | Rel residual | Sign correct? | Fits (<20%)? |
|---|---|---:|---:|---:|:---:|:---:|
| (a) c·ρ^n | c=0.00735, ρ=0.906 | +0.00408 | −0.00147 | 3.78 | **NO** | NO |
| (b) c·n^{-α}·(1/2)^n | c=5.74e-4, α=−3.87 | +0.00922 | −0.00147 | 7.27 | **NO** | NO |
| (c) c·log(n)·(1/2)^n | c=0.0429 | +0.001202 | −0.00147 | 1.82 | **NO** | NO |
| (d) c_1·ρ_1^n + c_2·ρ_2^n | a=1.008, b=−0.328, complex ρ, \|ρ\|=0.573 | +0.001611 | −0.00147 | 2.10 | **NO** | NO |
| (e) c·cos(ωn+φ)·(1/2)^n | c≈0.13, ω≈0.29, φ≈−1.4 | +0.001901 | −0.00147 | 2.29 | **NO** | NO |

**Every ansatz predicts the wrong sign at n=6.** The best (closest in magnitude) is ansatz (c) log, but it's still 182% off and on the wrong side. The threshold 20% is not approached by any of the five.

---

## Why all ansatze fail in the same direction

The n=2..5 fit window contains only positive δ values, monotone-decaying after the n=3 peak. Each ansatz form (single sign, smooth power-law/log/geometric extrapolation, or shallow-oscillation Prony recurrence) inherits this monotone-positive structure and predicts δ_6 in the positive ~+0.001 to +0.009 range. The actual δ_6 is negative (−0.00147). No ansatz captures the qualitative break.

This is a strong empirical signal that the analytic structure of δ_n cannot be reduced to any of the five pre-registered single-pattern candidates **on this data sample**. Either the sample is too small to capture the true oscillation period (the implied cosine in ansatz d / e has period ≈ 12.7 steps — we have 5 points), or δ_n has richer structure (three+ singularities, branch-cut + multiple poles, or accumulating singularities) not captured by any two-mode superposition.

Disposition follows in DISPOSITION.md.
