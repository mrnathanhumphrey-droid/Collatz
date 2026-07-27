# RESULT — P6H: the telescope collapses S to ONE limit — S_{i+1}=2·3ⁱ[4R_e(2)−R_e(0)] — and it's overshooting 7/30 (2026-07-26)

**Probe:** `probe_p6h.py`. Wilson's telescope: `T_i := 3ⁱ[4R_e⁽ⁱ⁾(2)−R_e⁽ⁱ⁾(0)] = Σ_k4⁻ᵏγ_i(k) = ½S^{(i)} = ½S_{i+1}`,
so `Λ_i = T_i−T_{i−1}` telescopes and the whole S-ladder is two autocorrelation values of the even sub-measure ν_e.

## (0) The identity is EXACT
`S_{i+1} = 2·3ⁱ[4R_e⁽ⁱ⁾(2)−R_e⁽ⁱ⁾(0)]` vs the certified S-ladder (`2/3 + 2Σ_{j≤i}Λ_j`, shellA): **diff 1e-15, i=1..6**.
`T_1 = 5/21`, `S_2 = 10/21` exact. The S-ladder has a closed form in two numbers per level from one sub-measure.
`S_∞ = 2·lim T_i`; `7/15 ⟺ lim T_i = 7/30 = 0.233333`.

## (1)(2) Two linearly-divergent collision sequences; the cancellation is NOT exact at depth
Both `3ⁱR_e(0)` and `4·3ⁱR_e(2)` diverge linearly (like `X_i = 3ⁱΣν²`, the ℓ² aggregate):
| quantity | slope (fit i=7..12) |
|---|---|
| `3ⁱR_e(0)` | 0.031211 |
| `4·3ⁱR_e(2)` | 0.031583 |
| **residual** | **+0.000373** |
| `X_i` | 0.156053 |
So `T_i = 4·3ⁱR_e(2) − 3ⁱR_e(0)` is the **residue of a near-cancellation** (Wilson's read, confirmed) — but the slopes
do **not** exactly cancel at this depth: a **residual +0.00037/level**, which *is* the current `Λ_i` (i=8..12 ≈ +0.00037).
The residual is **slowly shrinking** (0.000369 → 0.000336 over i=8→12) — consistent with an eventual, very deep turnover,
not resolved here. `4·slope(R_e(2)) = slope(R_e(0))` exactly is not a finite-depth identity; whether it holds
asymptotically **is** the 7/15 question restated.

## (3) T_i has OVERSHOT 7/30 and is rising — 7/15 needs a turnover
| i | T_i | Λ_i | T_i − 7/30 |
|---|---|---|---|
| 1 | 0.238095 | −0.095238 | +0.004762 |
| 2 | 0.230787 | −0.007308 | −0.002546 |
| 5 | 0.233084 | +0.000327 | −0.000249 |
| 8 | 0.233330 | +0.000369 | **−0.0000038** (crosses) |
| 10 | 0.234084 | +0.000391 | +0.000751 |
| 12 | 0.234807 | +0.000337 | **+0.001474** |

`T_i` is non-monotone early, then from i≈7 rises steadily, **crosses 7/30 at i≈8 and keeps climbing** to 0.23481 at
i=12 (+0.00147 above). `Σ_{i=2..12}Λ_i = −0.00329` — has moved **up** from the i=6 partial (−0.00535) because Λ₇..₁₂ are
positive — while the target is `−1/210 = −0.00476`. So at observable depth the sum is moving **away** from 7/15.

## Verdict — the identity is exact and beautiful; 7/15 is (still) conditional on a deep turnover
The constant is now a **single limit** of one explicit sequence built from `R_e(2), R_e(0)` — no sums, no tails, no free
maps. But `T_i` has overshot `7/30` and is rising with a persistent positive `Λ_i ≈ +0.00037` (the un-cancelled residual
slope). `7/15` requires that residual to vanish asymptotically and `T_i` to turn back down to `7/30` — **not observed at
i≤12** (extrapolating the rise gives `S_∞ ≈ 0.476`, above `7/15 = 0.4667`). This is the corpus's **log-periodic-bore /
conditional-7/15** finding (README R23, CROSSING/LOG-PERIODIC — turnover conjectured near r≈27–31) restated in the
clean two-autocorrelation-value language: the whole question is whether `3ⁱ[4R_e⁽ⁱ⁾(2)−R_e⁽ⁱ⁾(0)]` turns over. **Pen
(Wilson):** the asymptotic slope-cancellation `4·slope(R_e(2)) − slope(R_e(0)) → 0` and the turnover — the residual is
shrinking (0.000369→0.000336) but slowly. Not at stake: P6D/E/F/G identities (all exact), P1LVL, BRIDGE2, R1–R30.
