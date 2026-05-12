# PADE_NUMERICAL_DATA — Phase 1: assembled ε_n table k=1..13

**Date:** 2026-05-12. Wilson follow-up probe extending R77.6 / PADE_EXTENSION to numerical ε_7..ε_13.

## Source provenance

Exact rationals for k=1..6 cached in the project. Numerical floats for k=7..13 from prior power-iteration / scipy.eigs / FFT cross-checked sessions documented in `result_epsilon_*.md` and `probe_epsilon_12/`, `probe_epsilon_13/`.

| k | ε_k (float64) | method | source file |
|---|---|---|---|
| 1 | +2.0000000000e-01 = +1/5 | exact rational (cached) | result_77_6_pade_construction.py Stage 1 |
| 2 | +9.5238095238e-03 = +1/105 | exact rational (cached) | result_77_6_pade_construction.py Stage 1 |
| 3 | -5.0919863259e-03 = -5191/1019445 | exact rational (cached) | result_77_6_pade_construction.py Stage 1 |
| 4 | -2.4522582483e-03 = -11346676448406637/4627031617157687115 | exact rational (cached) | result_77_6_pade_construction.py Stage 1 |
| 5 | -1.1517469151e-03 | exact rational (cached, ~62-digit) | result_77_6_pade_construction.py Stage 1 |
| 6 | -4.9790566522e-04 | exact rational (cached, ~219-digit) | result_77_6_pade_construction.py Stage 1 |
| 7 | -1.1752368304e-03 | float64 power iter on K_7 (1458 states), scipy.eigs cross-check ag. 1e-15 | result_epsilon_7.md |
| 8 | -7.4554636729e-04 | float64 power iter on K_8 (4374), eigs+FFT cross-checks 3.33e-16..1.78e-15 | result_epsilon_8.md |
| 9 | -7.5202571564e-06 | float64 power iter on K_9 (13122), eigs+FFT cross-checks 9.44e-16 | result_epsilon_9.md |
| 10 | +7.2075091711e-04 | float64 matrix-free + eigs + FFT (39366 states); agreement 3.55e-15 | result_epsilon_10.md |
| 11 | +1.5019670121e-03 | float64 matrix-free + FFT (118098 states); FFT agreement 2.00e-15 | result_epsilon_11.md |
| 12 | +2.2747137206e-03 | float64 matrix-free + FFT (354294 states); FFT agreement 1.33e-15 | probe_epsilon_12/epsilon_12_findings.md |
| 13 | +2.9482473172e-03 | FFT on cached pi_13_truncated (v_max=60 trunc, error ≈ 8.7e-19) | probe_epsilon_13/epsilon_13_findings.md |

## Sign pattern + envelope

|n |sign(ε_n)| ε_n| |ε_n|·2^n |
|--|---|---|---|
|1 | + | 2.000e-01 | 4.0e-01 |
|2 | + | 9.524e-03 | 3.81e-02 |
|3 | − | 5.092e-03 | 4.07e-02 |
|4 | − | 2.452e-03 | 3.92e-02 |
|5 | − | 1.152e-03 | 3.68e-02 |
|6 | − | 4.979e-04 | 3.19e-02 |
|7 | − | 1.175e-03 | **1.504e-01** ← 4.7× jump |
|8 | − | 7.455e-04 | 1.91e-01 |
|9 | − | 7.520e-06 | 3.85e-03 ← near-zero node |
|10 | + | 7.208e-04 | 7.38e-01 |
|11 | + | 1.502e-03 | 3.08e+00 |
|12 | + | 2.275e-03 | 9.32e+00 |
|13 | + | 2.948e-03 | **2.42e+01** |

Sign sequence: `+ + − − − − − − − + + + +`. Single zero-crossing between k=9 and k=10. Envelope |ε_n|·2^n is NOT bounded — grows from 0.04 (n=5..6) to 24+ at n=13. The "n=2..6 plateau near 0.04" is the visible portion of a longer transient; n=7 enters a different regime (slow oscillating mode per STATE.md, ρ ≈ 0.984 in k-space, period ≈ 9.2).

## Coefficients for f̃(z) = (E(z) − ε_1 z)/z²

Working series c_j = ε_{j+2} for j = 0..11 (12 coefficients total, supporting Padé up to m+n ≤ 11):

| j | c_j = ε_{j+2} | float |
|---|---|---|
| 0 | ε_2 | +9.5238095238e-03 |
| 1 | ε_3 | -5.0919863259e-03 |
| 2 | ε_4 | -2.4522582483e-03 |
| 3 | ε_5 | -1.1517469151e-03 |
| 4 | ε_6 | -4.9790566522e-04 |
| 5 | ε_7 | -1.1752368304e-03 |
| 6 | ε_8 | -7.4554636729e-04 |
| 7 | ε_9 | -7.5202571564e-06 |
| 8 | ε_10 | +7.2075091711e-04 |
| 9 | ε_11 | +1.5019670121e-03 |
| 10 | ε_12 | +2.2747137206e-03 |
| 11 | ε_13 | +2.9482473172e-03 |

12 coefficients → Padé well-defined for m + n ≤ 11. Diagonal sequence [n/n] runs n=1..5 (since n=6 needs m+n=12, i.e. 13 coefficients of f̃ → ε_14 which we do not have). The probe builds diagonal up to [5/5] and near-diagonal [4/5], [5/4], [4/6], [6/4], [3/6], [6/3] etc as 11-budget allows.

**Correction:** Phase budget. With 12 coefficients we have m+n+1 ≤ 12, i.e. m+n ≤ 11. Diagonals up to [5/5] (uses 11), and [6/6] would need 13 coefficients (ε_14) — UNAVAILABLE. Reported diagonals: [1/1], [2/2], [3/3], [4/4], [5/5].
