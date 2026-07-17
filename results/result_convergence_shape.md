# Result: convergence shape of S_k -> 7/15 (q=3)

**Date:** 2026-05-05.  Characterizes the functional shape of S_k - 7/15 for k=1..5 using exact rational data cached in `result_q_sweep_test_1_envelope.csv`.

## Headline findings

1. **Sign pattern (corrected):** (+, +, -, -, -). One sign flip between k=2 and k=3; monotone-from-below for k>=3. The brief's premise of damped oscillation (model b) is **not supported** by the cached data — there is no oscillation, only a single sign flip. Brief's stated `S_5 ≈ 0.4673` and `eps_5 ≈ +0.0006` were wrong; correct values are `S_5 = 0.4655`, `eps_5 = -0.00115`.

2. **Envelope rate (log|eps_k| vs k):**
   - On k=1..5: rate ≈ 0.311 (k=1 transient dominates and inflates the rate).
   - **On k=2..5 only: rate = 0.493210, AICc = -11.1.** Essentially tight to 1/2 — the rate-1/2 envelope is confirmed in the asymptotic regime.

3. **Best signed fit: two-mode geometric with rates fixed at 1/2 and 1/3:**
   ```
       eps_k ≈ -0.4918 · (1/2)^k + 1.3268 · (1/3)^k
   ```
   AICc = **-35.977** (best of all candidates on linear-signed scale, by a large margin).
   The signs and rates have natural structural meaning:
   - `(1/2)^k` = 2-adic stripping rate (mean of `2^(-v)` under v ~ Geom(1/2) gives `1/3`, not `1/2`; the `1/2` here is the Markov chain second-eigenvalue scale)
   - `(1/3)^k` = 3-adic Bohr decay rate (matches `E[2^(-v)] = 1/3` from R66 framework)
   - The two contributions have *opposite signs* — exactly the structure that lets the sequence cross zero between k=2 and k=3.
   The amplitudes A ≈ -0.49 and B ≈ +1.33 are close to (but not exactly) -1/2 and +4/3; clean rationals would imply 7/36 at k=1 instead of the actual 1/5, so the closed form is approximate, not exact.

4. **Two-mode with FREE rates** is overparameterized (4 params, 5 points; AICc undefined). Fit degenerates to r1 ≈ r2 ≈ 0.045 (two redundant copies of one mode).

5. **Algebraic correction `A·r^k·k^α`** does not fit cleanly — best fit hits weird r > 1 with α ≈ -4.7. Not informative at this resolution.

## Bottom line

- **Rate**: confirmed 1/2 in asymptotic regime (k≥2). Pure geometric on log scale, AICc clean.
- **Shape**: best described by **two-mode (1/2)^k + (1/3)^k with opposite-sign amplitudes**. The sign flip between k=2 and k=3 falls out of this naturally — `(1/3)^k` dominates at small k (positive, large B), `(1/2)^k` takes over at larger k (negative, smaller |A|).
- **Note on AICc comparability**: envelope (log-scale) and signed (linear-scale) fits have different objectives; AICc values shown in the fits table aren't directly comparable across the two scales. Within each scale the best models are clear ((a-env-late) on log, (c-fixed) on linear).

## eps_k table

| k | S_k (decimal) | eps_k (rational) | eps_k (decimal) | sign | |eps_{k+1}/eps_k| |
|---|---|---|---|---|---|
| 1 | 0.6666666667 | 1/5 | +2.000000e-01 | + | 0.0476 |
| 2 | 0.4761904762 | 1/105 | +9.523810e-03 | + | 0.5347 |
| 3 | 0.4615746803 | -5191/1019445 | -5.091986e-03 | - | 0.4816 |
| 4 | 0.4642144084 | -11346676448406637/4627031617157687115 | -2.452258e-03 | - | 0.4697 |
| 5 | 0.4655149198 | -(59-digit num)/(62-digit den) | -1.151747e-03 | - |  |

## Functional-form fits

| model | AICc | SS_res | formula |
|---|---|---|---|
| (a-env) pure_geom_envelope | 5.717 | 2.1230e+00 | `|eps_k| = (2.5592e-01) * (0.3113)^k` |
| (a-env-late) pure_geom k>=2 | -11.147 | 4.5149e-03 | `|eps_k| = (4.0602e-02) * (0.4932)^k  (fit on k=2..5)` |
| (d-env) alg_correction | 15.451 | 2.7243e-01 | `|eps_k| = (9.3419e-02) * (1.9294)^k * k^(-4.6633)` |
| (c-fixed) two_mode 1/2 + 1/3 | -35.977 | 5.0752e-04 | `eps_k = (-4.9178e-01)*(1/2)^k + (+1.3268e+00)*(1/3)^k` |
| (c-free) two_mode free | inf | 3.7923e-05 | `eps_k = (+2.1065e+00)*(0.0450)^k + (+2.3320e+00)*(0.0451)^k  [OVERPARAMETERIZED: 4 params, 5 points]` |

**Pure-geom envelope rate (k=1..5):** 0.311273  (reference 1/2 = 0.500000)
**Pure-geom envelope rate (k=2..5 only, drops k=1 transient):** 0.493210

The k=2..5 fit is the asymptotic-regime estimate; the k=1 transient is a different scale and inflates the rate when included.

## Algebraic structure

- eps_1 = 1/5
- eps_2 = 1/105 = 1/(3·5·7)
- eps_3 = -5191/1019445 = -(29·179)/(3·5·7²·19·73)
- eps_2/eps_1 = 1/21 = 1/(3·7) — clean rational factor
- eps_3/eps_2 ≈ -0.535 — does not factor cleanly into small primes

Pattern observation: eps_1 and eps_2 have numerator 1 and denominators built from small odd primes. eps_3 onward has compound numerators with no clean structural factor. Suggests the small-k behavior (k = 1, 2) may have a closed-form origin while k >= 3 is genuinely algebraic-non-elementary at this computed precision.

## What was skipped

- Step 4(b) damped-osc fit: skipped (no oscillation to fit)
- Step 6 Hateley-style block-average comparison: skipped (no lift data in workspace)

## Files

- `result_convergence_shape.py` — script
- `result_convergence_shape_data.csv` — eps_k as exact rationals + decimals
- `result_convergence_shape_fits.csv` — model AICc / params table
- `result_convergence_shape.md` — this writeup