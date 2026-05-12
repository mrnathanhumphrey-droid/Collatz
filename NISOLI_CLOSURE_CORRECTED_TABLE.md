# NISOLI_CLOSURE_CORRECTED_TABLE — (r, A, K) tabulation at λ = 43/45

**Date:** 2026-05-12. Phase 2 of Route A Nisoli closure-inequality test at the corrected rate.

---

## 1. Parameters

- **M_3'' = 24.43** (exact, from Phase 1).
- **p = 3** (q = 3 in the c = 7/45 context).
- **Bilinear bound regimes:**
  - r ≤ 3: strict `|K_bil| = 2√N` (PATH2 family-level).
  - r ≥ 4: polylog-free `|K_bil| = 2√p · √N = 2√3 · √N` (HENSEL Hensel-digit-extraction).

The closure inequality: `|K_bil(r)| · K^{−A} · M_3'' < 1`.

---

## 2. Master table

Columns: r, N = 3^{r−1}, √N, |K_bil|, A (algebraic decay exponent), K (truncation level), K^{−A}, product = |K_bil|·K^{−A}·M_3'', fires?

| r | N | √N | \|K_bil\| | A | K | K^{−A} | product | fires |
|---|---|---|---|---|---|---|---|---|
| **3** | 9 | 3 | 6 | 1 | 6 | 0.1667 | 24.43 | **NO** |
| 3 | 9 | 3 | 6 | 1 | 10 | 0.1000 | 14.66 | NO |
| 3 | 9 | 3 | 6 | 1 | 20 | 0.0500 | 7.33 | NO |
| 3 | 9 | 3 | 6 | 1 | 50 | 0.0200 | 2.93 | NO |
| 3 | 9 | 3 | 6 | 1 | 100 | 0.0100 | 1.466 | NO |
| 3 | 9 | 3 | 6 | 1 | 150 | 0.00667 | 0.977 | **YES** |
| 3 | 9 | 3 | 6 | 2 | 6 | 0.02778 | 4.07 | NO |
| 3 | 9 | 3 | 6 | 2 | 10 | 0.01 | 1.466 | NO |
| 3 | 9 | 3 | 6 | 2 | 13 | 0.005917 | 0.867 | **YES** |
| 3 | 9 | 3 | 6 | 2 | 20 | 0.0025 | 0.366 | YES |
| 3 | 9 | 3 | 6 | 2 | 50 | 0.0004 | 0.0586 | YES |
| 3 | 9 | 3 | 6 | 2 | 100 | 0.0001 | 0.01466 | YES |
| 3 | 9 | 3 | 6 | 3 | 6 | 4.630e-3 | 0.679 | **YES** |
| 3 | 9 | 3 | 6 | 3 | 10 | 1.000e-3 | 0.1466 | YES |
| 3 | 9 | 3 | 6 | 3 | 20 | 1.250e-4 | 0.01833 | YES |
| 3 | 9 | 3 | 6 | 5 | 6 | 1.286e-4 | 0.01885 | YES |
| 3 | 9 | 3 | 6 | 10 | 6 | 1.654e-8 | 2.42e-6 | YES |
| **4** | 27 | 5.196 | 18 | 1 | 6 | 0.1667 | 73.3 | NO |
| 4 | 27 | 5.196 | 18 | 1 | 100 | 0.0100 | 4.40 | NO |
| 4 | 27 | 5.196 | 18 | 1 | 500 | 0.0020 | 0.880 | **YES** |
| 4 | 27 | 5.196 | 18 | 2 | 6 | 0.02778 | 12.21 | NO |
| 4 | 27 | 5.196 | 18 | 2 | 10 | 0.01 | 4.40 | NO |
| 4 | 27 | 5.196 | 18 | 2 | 21 | 0.002268 | 0.998 | **YES** |
| 4 | 27 | 5.196 | 18 | 2 | 50 | 0.0004 | 0.176 | YES |
| 4 | 27 | 5.196 | 18 | 3 | 8 | 1.953e-3 | 0.859 | **YES** |
| 4 | 27 | 5.196 | 18 | 3 | 10 | 1.000e-3 | 0.440 | YES |
| 4 | 27 | 5.196 | 18 | 3 | 20 | 1.250e-4 | 0.0550 | YES |
| 4 | 27 | 5.196 | 18 | 5 | 6 | 1.286e-4 | 0.0566 | YES |
| 4 | 27 | 5.196 | 18 | 10 | 6 | 1.654e-8 | 7.27e-6 | YES |
| **5** | 81 | 9 | 31.18 | 1 | 6 | 0.1667 | 126.9 | NO |
| 5 | 81 | 9 | 31.18 | 1 | 1000 | 0.001 | 0.762 | **YES** |
| 5 | 81 | 9 | 31.18 | 2 | 6 | 0.02778 | 21.16 | NO |
| 5 | 81 | 9 | 31.18 | 2 | 28 | 1.276e-3 | 0.972 | **YES** |
| 5 | 81 | 9 | 31.18 | 2 | 50 | 0.0004 | 0.305 | YES |
| 5 | 81 | 9 | 31.18 | 3 | 10 | 1.000e-3 | 0.762 | **YES** |
| 5 | 81 | 9 | 31.18 | 3 | 20 | 1.250e-4 | 0.0953 | YES |
| 5 | 81 | 9 | 31.18 | 5 | 6 | 1.286e-4 | 0.0980 | YES |
| 5 | 81 | 9 | 31.18 | 5 | 10 | 1.000e-5 | 7.62e-3 | YES |
| **6** | 243 | 15.59 | 54.0 | 1 | 6 | 0.1667 | 219.9 | NO |
| 6 | 243 | 15.59 | 54.0 | 2 | 6 | 0.02778 | 36.65 | NO |
| 6 | 243 | 15.59 | 54.0 | 2 | 37 | 7.305e-4 | 0.964 | **YES** |
| 6 | 243 | 15.59 | 54.0 | 3 | 11 | 7.513e-4 | 0.992 | **YES** |
| 6 | 243 | 15.59 | 54.0 | 5 | 6 | 1.286e-4 | 0.1697 | YES |
| 6 | 243 | 15.59 | 54.0 | 5 | 7 | 5.949e-5 | 0.0785 | YES |
| 6 | 243 | 15.59 | 54.0 | 10 | 6 | 1.654e-8 | 2.18e-5 | YES |

---

## 3. Minimum-firing K at each (r, A)

`K^A > |K_bil(r)| · M_3''` ⟹ `K_min = (|K_bil|·M_3'')^{1/A}`. Minimum-firing K at each cell (assuming `C_A = 1` in the Tao Prop 1.17 statement):

| r | \|K_bil\| | thr = \|K_bil\|·M_3'' | A=1 | A=2 | A=3 | A=5 | A=10 |
|---|---|---|---|---|---|---|---|
| 2 | 3.464 | 84.6 | 85 | 9.20 | 4.39 | 2.43 | 1.56 |
| **3** | 6 | **146.6** | 147 | **12.11** | **5.27** | 2.70 | 1.64 |
| 4 | 18 | 439.7 | 440 | 20.97 | 7.61 | 3.36 | 1.83 |
| 5 | 31.18 | 761.7 | 762 | 27.60 | 9.14 | 3.78 | 1.94 |
| 6 | 54.0 | 1,319 | 1,319 | 36.32 | 10.97 | 4.20 | 2.05 |
| 8 | 162.0 | 3,957 | 3,957 | 62.91 | 15.81 | 5.20 | 2.28 |
| 10 | 486.0 | 11,873 | 11,873 | 108.96 | 22.80 | 6.46 | 2.54 |

Reading: at r = 3 (the strict-2√N regime), closure fires already at A = 2 with K ≥ 13, or A = 3 with K ≥ 6.

---

## 4. Cross-comparison to rate-1/2 (M3_CLOSURE_TABLE.md)

`M3_CLOSURE_TABLE.md` tabulated the **same** inequality at rate-1/2, where the falsified T_3 had M_3 = 100 (anticipated) or 1000 (R77.2 loose):

| Quantity | rate-1/2 (M3_CLOSURE_TABLE M_3=100) | rate-1/2 (M_3=1000) | corrected rate 43/45 (M_3''=24.43) |
|---|---|---|---|
| min A at K=6, r=3 | 2.65 | 3.94 | 3 (closes at K=6 with A=3 product 0.679) |
| min A at K=10, r=3 | 2.06 | 3.06 | ~2.18 (closes at A=3 product 0.147) |
| min A at K=20, r=3 | 1.58 | 2.36 | ~1.80 |

The corrected-rate closure is **EASIER than the rate-1/2 loose case (M_3=1000) and slightly easier than the rate-1/2 tight case (M_3=100)** for r=3 because M_3'' = 24.43 < 100 < 1000.

This is the counterintuitive finding from Phase 1: even though λ = 43/45 is closer to 1, T_lead's 2×2 rank-1 structure (κ ≈ 1.086) makes M_3'' much smaller than the rate-1/2 T_3's M_3 (which had κ ≈ 50-1000 due to companion-matrix non-normality).

---

## 5. Pre-registered hypotheses — Phase 2 update

| Hypothesis | Tabulation says |
|---|---|
| H_NISOLI_CLOSES_AT_R3 (firing at r=3 with reasonable A) | **CANDIDATE** — closure fires at (r=3, A=3, K≥6) with product 0.679, or (r=3, A=2, K≥13). |
| H_NISOLI_CLOSES_AT_LARGER_R (firing at r=4 or 5) | **CANDIDATE** — fires at (r=4, A=3, K≥8), (r=4, A=2, K≥21), (r=5, A=3, K≥10). |
| H_NISOLI_NEEDS_TIGHTENING | RULED OUT for the M_3'' bound itself — already tight at 24.43. |
| H_NISOLI_DOESNT_CLOSE | RULED OUT — closure clearly fires at many cells. |
| H_A_EXTRACTION_HARD | **PENDING** — Phase 4. Cells fire if A is givable; Tao's A is **free parameter** but `C_A` is INFEASIBLE per BOOKKEEPING_PHASE1_DISPOSITION. |

Phase 3 (NISOLI_CLOSURE_CORRECTED_FIRING) documents the smallest firing cells. Phase 4 (HONEST_READ) addresses whether the A those cells require is actually deliverable from the project's machinery.

---

## 6. Tao Prop 1.17 Pareto frontier (A vs implicit C_A)

The Nisoli closure inequality assumes the algebraic-decay form `|μ̂_n(ξ)| ≤ K^{−A}` (i.e., C_A = 1). In reality, Tao Prop 1.17 says `|μ̂_n(ξ)| ≤ C_A · n^{−A}`, where:

- Larger A in Tao's proof requires larger `C_A` (the constant absorbs the bookkeeping cost).
- BOOKKEEPING_PHASE1_DISPOSITION §"Highest-difficulty constants" shows `C_A` grows at least like `exp(exp(A²))` under faithful bookkeeping, and at the optimistic floor `A^{O(A)}` under any reading.

So the correct closure inequality is:

> `|K_bil| · C_A · K^{−A} · M_3'' < 1`

with `C_A = exp(exp(A²))` or `A^{O(A)}` floor. The C_A factor dominates: for closure to fire we need
- `K^A > |K_bil| · C_A · M_3''`

At r = 3, A = 3, K = 6: required `6^3 = 216 > 146.6 · C_3`. So `C_3 < 1.47`. Under any Tao bookkeeping bound, `C_3 ≥ 3^{O(3)} = 27` or worse. **The optimistic floor already blows the budget.**

This is the same finding as M3_CLOSURE_TABLE §"caveats A2" and BOOKKEEPING_PHASE1_DISPOSITION §"Highest-difficulty constants" §C_A failure: **the conditional `C_A = O(1)` parameterisation rules in the tabulation are NOT what Tao's proof delivers.**

---

## 7. Verdict for Phase 2

**The closure inequality fires at many (r, A, K) cells** under the optimistic `C_A = 1` parameterisation. The smallest firing cells are:

- **r=3, A=3, K=6**: product 0.679 (closure fires; M_3'' · |K_bil| · 6^{−3} = 24.43 · 6 · 1/216 = 0.679)
- r=3, A=2, K=13: product 0.867
- r=4, A=3, K=8: product 0.859
- r=4, A=2, K=21: product 0.998

**HOWEVER**, the realistic Tao Prop 1.17 — with the actual super-exponential `C_A` — does **not** deliver `A ≥ 2` with `C_A · K^{−A} < (|K_bil|·M_3'')^{−1}` at any K ≤ 10. Phase 4 documents this honestly.

The Phase 2-3 tabulation gives **conditional firing** ("closure fires IF Tao delivers polynomial-in-A constants, which it does not"). The unconditional disposition routes to Phase 4.

---

## 8. Files

- nisoli_closure_corrected.py — main-thread verification script
- NISOLI_CLOSURE_CORRECTED_PHASE1.md — closure inequality articulated
- NISOLI_CLOSURE_CORRECTED_TABLE.md (this file)
- NISOLI_CLOSURE_CORRECTED_FIRING.md — Phase 3 firing-cells documentation
- NISOLI_CLOSURE_CORRECTED_HONEST_READ.md — Phase 4 A-extractability
