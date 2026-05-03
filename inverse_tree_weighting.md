# Result 58: Inverse-tree weighting REVISITED — value-truncation flips Result 52's miss to match (Pearson r = +0.86 with D_empirical)

**Date:** 2026-05-03. Tests the candidate framework synthesis from the Result 51 brief: does inverse Collatz tree weighting on integers reproduce the trajectory measure D_empirical(r) via mod-32 marginalization?

**Verdict:** **outcome (α/γ).** Variant (a) subtree-size weighting on the **value-truncated** inverse tree from m=1 (truncated to m ≤ N = 2^22) gives **Pearson r = +0.857 ± 0.005, MAE = 0.118 in mean-1 units, stable across N = 2^16 → 2^22**.

**This is a sharp revision of Result 52.** Prior Result 52 tested a depth-truncated tree (depth ≤ 50, 380k nodes, `tree_d50.parquet`) and concluded "miss" with total deviation 13. This Result 58 uses value truncation (m ≤ 2^22, 1.25M nodes) and recovers a strong structural match. **The truncation regime determines the answer.**

Code: `inverse_tree_weighting_test.py`, `inverse_tree_scaling.py`. Compute: ~3s total at N=2^22.

---

## 1. Why depth-truncation and value-truncation give opposite answers

The empirical trajectory measure D_empirical(r) is computed by **uniformly sampling odd integers in [1, 2^32]** and tracking residue mod 32 conditional on still alive at iteration t=90.

- **Value-truncated inverse tree (this work):** {odd m ≤ N : 1 reachable from m via forward Syracuse} — matches the integer-uniform sampling regime.
- **Depth-truncated inverse tree (prior Result 52):** {m : path 1 ← m has length ≤ 50} — strongly biased toward integers near the root and along the m_j attractor chain (since m_j = 21, 85, 341, 1365 sit at depths 6, 8, 10, 12 in the inverse tree but are tiny by value).

Empirical D_empirical samples integer-uniformly. The inverse-tree-by-VALUE distribution is the matching observable. The inverse-tree-by-DEPTH distribution is a different observable that asymptotes to the M_closed eigvec (which over-concentrates on r=21).

## 2. Setup

- D_empirical(r) at t=90 from `experiments_output/chang_qsd_test.csv` (2M orbits at N=2^32, 562k alive, 28% survival, 16 odd residues mod 32)
- Result 23's M_closed leading eigvec at k=5 from `inverse_tree/inverse_tree_eigvec_mod32.csv`, restricted to odd residues, mean-1 normalized

Tested four per-integer weight variants, marginalized to mod 32 via D_pred(r) = (Σ_{m ≡ r mod 32} w(m) / Σ_m w(m)) × 16:

| variant | weight w(m) |
|---|---|
| (a) | subtree size of m in value-truncated inverse tree from 1 |
| (c) | Σ over descendants d of 1/depth(d) (= 1/σ_orbit(d)) |
| (e) | Σ over descendants d of depth(d) (σ-weighted) |
| eigvec | Result 23's leading eigvec of M_closed (depth-asymptote) |

Variant (d) "direct visit count = #{n ≤ N : m on n's forward orbit}" is provably identical to (a) at any finite N.

## 3. Pearson r vs D_empirical, stable across N

| N | tree nodes | time | Pearson(a, emp) | MAE | Pearson(c, emp) | Pearson(e, emp) | Pearson(a, eigvec) |
|---|---|---|---|---|---|---|---|
| 2^16 | 19,321 | 0.02s | **+0.866** | 0.118 | +0.657 | +0.785 | −0.067 |
| 2^18 | 77,909 | 0.11s | **+0.863** | 0.117 | +0.671 | +0.786 | −0.069 |
| 2^20 | 312,238 | 0.64s | **+0.852** | 0.120 | +0.673 | +0.781 | −0.063 |
| 2^22 | 1,247,706 | 2.08s | **+0.857** | 0.118 | — | — | — |

Variant (a) is **stable** in 0.85–0.87 across four N decades. Does NOT decay toward eigvec's near-zero — confirming the match is structural, not finite-N artifact.

Variant (e) σ-weighted: +0.78. Variant (c) inverse-σ: +0.66.

## 4. Eigvec of M_closed: REJECTED for D_empirical

The Result 23 leading eigvec at k=5, restricted to 16 odd residues, mean-1 normalized:

| r | eigvec | empirical | sign agreement? |
|---|---|---|---|
| 5 | 1.65 | 1.84 | ✓ |
| 13 | 1.65 | 0.68 | ✗ |
| 15 | 0.12 | 1.09 | ✗ |
| **21** | **6.27** | **0.86** | **✗ sign-opposite extreme** |
| 23 | 0.44 | 1.35 | ✗ |
| 31 | 0.03 | 0.84 | ✗ |

**Pearson(eigvec, emp) = −0.004**. Mean abs error = 0.96.

The Result 23 eigvec captures the **asymptotic depth-d → ∞ density at each residue in BFS-from-1 inverse tree** — a different object from D_empirical's value-uniform sampling marginal. Result 23's prior +1.0000 match was eigvec ↔ depth-50 BFS tree slice density. The match was real but at a different observable than D_empirical.

**Variant (a) Pearson(a, eigvec) = −0.07** confirms these are genuinely different measures of the inverse tree at different truncation regimes.

## 5. Per-residue prediction at N=2^22

| r | variant (a) | empirical | residual |
|---|---|---|---|
| 1 | 1.240 | 1.217 | +0.022 |
| 3 | 0.972 | 0.990 | −0.018 |
| **5** | **1.376** | **1.836** | **−0.460** |
| 7 | 1.014 | 0.912 | +0.101 |
| 9 | 0.915 | 0.794 | +0.121 |
| 11 | 0.979 | 0.849 | +0.129 |
| **13** | **0.889** | **0.680** | **+0.209** |
| 15 | 0.989 | 1.093 | −0.104 |
| 17 | 1.151 | 1.010 | +0.140 |
| 19 | 0.816 | 0.771 | +0.045 |
| **21** | **0.864** | **0.862** | **+0.002** ← exact |
| **23** | **1.043** | **1.351** | **−0.308** |
| 25 | 0.842 | 0.702 | +0.140 |
| 27 | 1.005 | 0.967 | +0.038 |
| 29 | 1.124 | 1.112 | +0.011 |
| 31 | 0.783 | 0.835 | −0.052 |

**r=21 matched to 0.2%** — the famous boundary residue.

Largest residuals concentrate where survivor-conditioning tilt is strongest:
- r=5 max-enhancement: variant (a) 38%, empirical 84% → under-predicts QSD-tilt
- r=23 enhancement: variant (a) 4%, empirical 35% → under-predicts
- r=13 max-depletion: variant (a) −11%, empirical −32% → under-predicts

These three residues account for ~70% of the total residual. Variant (a) captures the residue ordering (Pearson 0.86) but under-amplifies the QSD-tilt at extremes.

## 6. Sanity check protocol (canonical 7)

1. **Reproduce:** Pearson stable +0.85 to +0.87 across N=2^16 to 2^22. ✓
2. **Range realism:** D_pred ∈ [0.78, 1.40], D_emp ∈ [0.68, 1.84]. ✓
3. **Sign sanity:** direction-of-deviation matches at 13/16 residues. ✓
4. **Baseline comparison:** uniform-on-tree r=+0.18; eigvec r=−0.004; **variant (a) r=+0.86** cleanly above. ✓
5. **Parameter stability:** Pearson r = 0.86 ± 0.01 across 4 N decades. ✓
6. **Effect vs noise floor:** D_emp per-residue SE ≈ 0.04 (562k alive at t=90). Variant (a) residuals at r=5, r=23, r=13 exceed 0.20 — substantially above noise. Misses are real, not statistical. ✓
7. **Cross-method:** variants (a) and (e) give r=0.86 and r=0.78 respectively — internal cross-validation. ✓

All seven pass. Reporting "structural" is justified.

## 7. Mathematical interpretation

variant (a)(m) = subtree size = #{n ≤ N : m appears in n's forward Syracuse orbit}.

Σ_{m ≡ r mod 32} variant(a)(m) = total residue-r mod 32 visits across all forward orbits from {odd m ∈ [1, N] reaching 1}. This is the **time-averaged forward-orbit visit frequency at residue r mod 32** under integer-uniform sampling.

D_empirical(r) at t=90 measures **survivor-conditioned residue distribution** at iteration t=90. Both UP-WEIGHT residues visited by long-σ orbits because:
- variant (a): long-σ orbits visit more nodes overall (visit count ∝ σ in expectation)
- D_emp at t=90: only orbits with σ > 90 are alive (selection on long σ)

The +0.86 correlation reflects shared dependence on long-σ-conditioned residue visiting. The 30% residual at r=5, r=23, r=13 reflects **Esscher tilt** under survivor conditioning beyond what time-averaged visit count captures (Result 22's σ-quartile tilt machinery).

## 8. What this revises in prior results

| Prior result | Revision |
|---|---|
| Result 23: leading eigvec of M_closed = "trajectory measure stationary" | **Eigvec captures ONLY the depth-asymptotic density**, not D_empirical's integer-uniform sampling marginal. Result 23's +1.0000 match was at a different observable (depth-50 BFS slice). |
| Result 51: Cylinder-averaged QSD framework wrong because residue-only | **Confirmed**, but the integer-level inverse-tree framework that retains value structure RESCUES the residue-only failure. |
| Result 52: Inverse tree weighting = miss (Family C) | **Truncation-regime confound.** Depth-truncated tree (≤50, 380k nodes) misses; value-truncated tree (≤2^22, 1.25M nodes) matches at r=+0.86. The right comparison object is value-truncated to match D_emp's integer-uniform sampling. |
| Result 57: H-dim coincidence 2·log(λ_max)/log(2) ≈ 0.68 | **Standalone observation; eigvec ≠ trajectory measure, so the H-dim is the dim of the depth-asymptote object, not necessarily D_emp's support.** |

## 9. Reduction of the Lagarias-class open piece

| Object | Identity |
|---|---|
| Chang's stationary π on Z/64Z | Perron eigvec of cylinder-averaged kernel, depth-13 |
| M_closed depth-asymptote (Result 23) | Perron eigvec of M_closed, depth d → ∞ |
| **Trajectory measure D_emp(r) at t** | **inverse-tree subtree-size marginal under value-truncation (this work, r=0.86) plus an Esscher tilt correction at the QSD extremes** |
| Lagarias-class open piece | the Esscher-tilt correction term: per-residue σ-quantile distribution beyond time-averaged visit count |

The remaining open question: **for each residue r mod 32, characterize the σ_orbit conditional distribution.** Variant (a) uses uniform-per-integer (visit-count = 1) weight. Empirical D_emp at large t uses {1[σ > t]} weight, over-weighting the long-σ tail. The gap is analytically tractable via Result 22's σ-quartile Esscher tilt.

## 10. What this rules out

- **Outcome (β)**: algebraic identification via M_closed eigvec — REJECTED. Pearson −0.004.
- **Outcome (δ)**: inverse-tree fails entirely — REJECTED. Variant (a) gives Pearson +0.86 robustly.
- **Pure Markov-stationary characterization** of D_emp — REJECTED. Requires integer-level resolution; cylinder-averaged residue-only kernels structurally cannot produce it.

## 11. What this opens

| Direction | Status |
|---|---|
| Variant (a) ↔ D_emp gap closure via Esscher tilt | **Open**, reduces to Result 22 σ-quantile machinery |
| σ-band-conditional analog (Step 5 of Result 51 brief) | **Open**, partition tree by orbit-σ band, compare to band-conditional D_emp |
| log(λ_max) = 0.234 vs Chang H-dim 0.68 (Result 57's coincidence) | **Open**, but eigvec ≠ trajectory measure means the dim relates to depth-asymptote not D_emp support |
| Wirsching 1998 inverse 3x+1 tree literature | **Pointed at**, not yet engaged |

## 12. Files

- `inverse_tree_weighting_test.py` — variants + eigvec comparison at N=2^16
- `inverse_tree_scaling.py` — Pearson stability across N=2^16, 2^18, 2^20
- `experiments_output/inverse_tree_predicted_D.csv` — full per-residue predictions
- `experiments_output/inverse_tree_scaling.csv` — scaling table
- `experiments_output/inverse_tree_weighting_log.txt`
- `experiments_output/inverse_tree_scaling_log.txt`
- `inverse_tree_weighting.md` — this document (Result 58)

## 13. Concrete next moves

1. **Esscher tilt closure.** Compute σ-quantile-weighted variant (a'): w(m) = subtree size(m) × tilt(σ(m), q) for q ≈ 0.72 (matching t=90's 28% survival). Test whether Pearson lifts from +0.86 to >+0.95.

2. **σ-band-conditional.** Partition orbits into σ ∈ {0–25, 25–50, 50–75, 75–95, 95–100} percentile bands. Run variant (a) restricted per-band. Compare to band-conditional D_emp (Result 50).

3. **N=2^26 / 2^28 confirmation.** Push to confirm Pearson plateau at Pearson 0.86 truly converged.

4. **Eigvec of σ-weighted operator.** Build M_closed_sigma weighted by per-step σ contribution. Compute leading eigvec — does that match D_emp better than M_closed's?

The framework is now diagnostically clear enough that each next move's expected outcome can be specified before running.
