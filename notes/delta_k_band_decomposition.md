# ΔK_band decomposition — methodology cleanup + Result 32 baseline correction (Result 33)

**Status.** Outcome between (a) and (b). The U-shape ΔK_band reported in
Result 25 was largely a methodology + E_V baseline artifact. After Result 32
cleanup (using per-orbit E_V_orbit instead of per-step E_per_step in K_bulk),
the boundary gap shrinks by an order of magnitude. q=0.875 closes within
bootstrap (0.3σ); other bands have residuals ±0.7, no longer monotone U-shape.

## The two K_eff measurements differ systematically

Within-band OLS (slope of σ on log_n given band) and first-passage
(slope of mean R on mean log threshold across 4 thresholds, restricted to
band) measure different objects on the same orbits.

5-seed bootstrap at N=2³⁶, 500k orbits:

| q     | K_ols mean ± sd | K_fp mean ± sd | gap (K_ols − K_fp) | sigma |
|------:|----------------:|---------------:|-------------------:|------:|
| 0.125 |  9.107 ± 0.151 |  6.929 ± 0.018 | +2.178 ± 0.152 | 14× |
| 0.375 | 10.451 ± 0.104 |  8.991 ± 0.029 | +1.460 ± 0.104 | 14× |
| 0.625 | 10.496 ± 0.104 | 10.624 ± 0.022 | −0.127 ± 0.121 |  1× |
| 0.875 | 12.030 ± 0.175 | 14.626 ± 0.071 | −2.595 ± 0.213 | 12× |
| 0.975 | 11.959 ± 1.275 | 18.891 ± 0.078 | −6.933 ± 1.306 |  5× |

10–14σ systematic differences. The U-shape from K_ols−K_fp matches the
brief's "ΔK_band U-shape" magnitude — confirming the U-shape lives largely
in the methodology choice.

K_fp(q=0.875) = 14.63 reproduces Result 23's plateau exactly. K_fp is the
"correct" first-passage K_eff from the prior threshold framework.

## The σ-identity decomposition is exact for K_ols

Algebraic identity: K_eff_ols(band) = (1+E_V)·slope_T + E_T·slope_V + κ_TVZ/Var(log_n).

Exact to machine precision (gap ≈ 0) across all bands and N — this is
σ = T·(1+V_orbit) tautology applied to OLS slope decomposition.

For middle bands, slope_T ≈ K_h/(1+E_V) closes within 3% via per-band
Jensen approximation under tight V|band:

| q     | slope_T | K_h/(1+E_V) | gap (×(1+E_V)) |
|------:|--------:|------------:|---------------:|
| 0.125 |  2.968 |  3.133 | −0.550 |
| 0.375 |  3.487 |  3.396 | +0.281 |
| 0.625 |  3.504 |  3.516 | −0.036 |
| 0.875 |  4.096 |  3.636 | +1.320 |
| 0.975 |  4.063 |  3.714 | +0.980 |

q=0.625 closes (−0.036). Other bands have K_h/(1+E_V) deviation of
±0.3–1.3, dominant in tails.

## The brief's ΔK_band U-shape was a Result 32 artifact

Brief defined ΔK_band = K_bulk(E_band) − K_eff_band, with K_band measured by
first-passage and E_band = per-step Esscher mean (Result 25).

But Result 32 showed E_per_step ≠ E_V_orbit — they differ by
−Cov[T,V|band]/E[T|band], which is structurally significant in the q=0.125
tail (−0.04) and the q=0.875–0.975 upper tails.

K_bulk(v) = (1+v)/(v·log2 − log3) is non-linear in v. Substituting per-step
v overestimates K_bulk in upper bands where E_per_step < E_V_orbit.

Comparison at N=2³⁶:

| q     | K_fp  | K_bulk(E_per_step) | ΔK_old | K_bulk(E_V_orbit) | ΔK_new |
|------:|------:|-------------------:|-------:|------------------:|-------:|
| 0.125 |  6.93 |  7.30 | +0.37 |  6.46 | **−0.47** |
| 0.375 |  8.99 |  9.82 | +0.83 |  9.11 | **+0.12** |
| 0.625 | 10.62 | 11.96 | +1.33 | 11.24 | **+0.61** |
| 0.875 | 14.63 | 17.85 | +3.22 | 14.60 | **−0.02** |
| 0.975 | 18.89 | 27.85 | +8.96 | 18.17 | **−0.72** |

**Reduction factor: up to 14× (q=0.975, +8.96 → −0.72).** q=0.875 closes
within bootstrap (−0.02 ± 0.07, 0.3σ). Middle band q=0.375 within 0.13.
The original "U-shape monotone-rising-to-+8.96" pattern disappears entirely.

## Closure status for constant 4

**Bulk.** Esscher per-step + Result 32 algebraic Cov[T,V|band]/E[T|band]
correction. Closed.

**Boundary (after R32 cleanup).** ΔK_band(q) = K_bulk(E_V_orbit(q)) − K_fp(q):
- q=0.125: −0.47
- q=0.375: +0.12
- q=0.625: +0.61
- q=0.875: −0.02 (closes within bootstrap)
- q=0.975: −0.72

±0.7 residual, no monotone trend, q=0.875 already closes. Outcome (a)-adjacent.

The remaining ±0.7 residual likely closes via either:
- Further Esscher refinement (per-band T-distribution beyond stationary mean)
- N-finite boundary effects (residual decays with N — testable)

NOT Lagarias-class — does not unify with constant 3 (per-j W_j) wall.

## What this means for v3.6

The "constant 4 boundary U-shape" framing from Result 25 substantially
overstates the open piece. With the correct E_V baseline:
- 4 of 5 bands within ±0.72
- 1 band closes within bootstrap (q=0.875)
- No monotone trend across bands

The substantive open piece is the ±0.7 residual structure, which is much
smaller than the original "+0.20 to +8.96" range and may close via finite-N
analysis or refined Esscher.

v3.6 framing: constant 4 closes structurally pending small (±0.7)
boundary-residual derivation; constant 3 (per-j W_j ⟺ ⟨σ_S | j⟩) remains
the single Lagarias-class open piece across the bridge constants.

## Files

- `experiments/65_delta_k_band_decomp.py` — σ-identity decomposition of K_ols
- `experiments_output/65_delta_k_band_decomp.csv`
- `experiments_output/65_delta_k_band_decomp_log.txt`
- `experiments/66_K_eff_method_compare.py` — K_ols vs K_fp on same orbits
- `experiments_output/66_K_eff_method_compare.csv`
- `experiments_output/66_K_eff_method_compare_log.txt`
