# D2 — BMT / bigraph predictions vs Syracuse numerical values

**Date:** 2026-05-15
**Mode:** E.

---

## 1. Syracuse measured moments (Reading B, level n=3, (Z/27)*, V_TRUNC=16)

Source files: `C:/Collatz/experiments_output/monotone_diagnostic_n3.json` (row d) and `C:/Collatz/experiments_output/n4_alternating_diagnostic.json` (rows b, f).

| Row | Moment | tr_pi | vac_pi | delta_1 | sum_entries |
|---|---|---|---|---|---|
| (b) | ϕ(X̃_1 · X̃_2) | 1.070×10⁻⁵ | 5.140×10⁻⁵ | 8.571×10⁻⁶ | 1.076×10⁻⁷ |
| (d) | ϕ(X̃_1 · X̃_2 · X̃_1) | 5.085×10⁻⁴ | 2.483×10⁻⁴ | 2.406×10⁻³ | **1.078×10⁻¹** |
| (f) | ϕ(X̃_1 · X̃_2 · X̃_1 · X̃_2) | **5.357×10⁻²** | **4.775×10⁻³** | **5.742×10⁻²** | **6.089×10⁻¹** |

(b) values are all at the noise floor (10⁻⁵ to 10⁻⁷). (d) and (f) are structurally non-zero per `D1_DISPOSITION.md`.

---

## 2. BMT predictions (every digraph)

| Row | Boolean | Monotone | Anti-mono | Tensor |
|---|---|---|---|---|
| (b) | 0 | 0 | 0 | 0 |
| (d) | 0 | 0 | 0 | 0 |
| (f) | 0 | 0 | 0 | ϕ(X̃_1²)·ϕ(X̃_2²) |

**See `D2_BMT_BIGRAPH_PREDICTIONS.md` §1 for derivations.**

## 3. Bigraph predictions (every pairwise relation)

| Row | Boolean | Monotone | Anti-mono | Free | Tensor |
|---|---|---|---|---|---|
| (b) | 0 | 0 | 0 | 0 | 0 |
| (d) | 0 | 0 | 0 | 0 | 0 |
| (f) | 0 | 0 | 0 | 0 | ϕ(X̃_1²)·ϕ(X̃_2²) |

**See `D2_BMT_BIGRAPH_PREDICTIONS.md` §2 for derivations.**

## 4. Row-by-row comparison

### Row (b) — Syracuse ≈ 0; BMT/bigraph predict 0 (every choice).

| Reduction | Syracuse | BMT prediction | Bigraph prediction | Match? |
|---|---|---|---|---|
| tr_pi | 1.070×10⁻⁵ | 0 | 0 | ✓ noise-floor agreement |
| vac_pi | 5.140×10⁻⁵ | 0 | 0 | ✓ |
| delta_1 | 8.571×10⁻⁶ | 0 | 0 | ✓ |
| sum_entries | 1.076×10⁻⁷ | 0 | 0 | ✓ |

**Verdict at row (b): every BMT/bigraph choice matches Syracuse.** ✓

### Row (d) — Syracuse measurably non-zero; BMT/bigraph predict 0 (every choice).

| Reduction | Syracuse | BMT prediction | Bigraph prediction | Match? |
|---|---|---|---|---|
| tr_pi | 5.085×10⁻⁴ | 0 | 0 | **✗** (≥ 30× control noise 1.7×10⁻⁶) |
| vac_pi | 2.483×10⁻⁴ | 0 | 0 | **✗** (≥ 200× control noise) |
| delta_1 | 2.406×10⁻³ | 0 | 0 | **✗** (≥ 400× control noise) |
| sum_entries | **1.078×10⁻¹** | 0 | 0 | **✗** (≥ 7500× control noise) |

**Verdict at row (d): NO BMT digraph and NO bigraph configuration predicts the measured Syracuse non-vanishing.** ✗

The row (d) mismatch is the **decisive gap**: BMT and bigraph both predict 0 because the peak / interval position (position 2 with label j_2) is ALWAYS isolated in kerG (BMT) or as a singleton block (bigraph). The centered factor ϕ(X̃_{j_2}) = 0 kills the prediction. But Syracuse's measured value 0.108 (sum_entries) corresponds to a structural cross-coupling through the b_prior-dependent phase of X̃_{j_2}, which is a Reading-B artifact NOT captured by either framework's moment formula.

### Row (f) — Syracuse non-zero; only TENSOR predicts non-zero.

| Reduction | Syracuse | BMT(mono) | BMT(tensor) | Bigraph(tensor) |
|---|---|---|---|---|
| tr_pi | 5.357×10⁻² | 0 ✗ | ϕ(X̃_1²)·ϕ(X̃_2²) (need numeric) | same |
| vac_pi | 4.775×10⁻³ | 0 ✗ | ϕ(X̃_1²)·ϕ(X̃_2²) | same |
| delta_1 | 5.742×10⁻² | 0 ✗ | ϕ(X̃_1²)·ϕ(X̃_2²) | same |
| sum_entries | **6.089×10⁻¹** | 0 ✗ | ϕ(X̃_1²)·ϕ(X̃_2²) | same |

The auditor predicted BMT under "monotone digraph" gives ker_G = {{1,3},{2,4}}. **This is incorrect** under the verbatim BMT Defn 2.8 + Prop 3.5(iii). The verbatim "monotone digraph" E = {(i,j) : j < i} has edge (j_2, j_1) (from larger to smaller) but NOT (j_1, j_2). The 2 ∼ 4 kernel-block requires (i_3, i_2) = (j_1, j_2) ∈ E — which is the WRONG direction for monotone. So position 4 stays a singleton, and BMT-monotone gives 0.

**Only the TENSOR case (both directions in E) gives kerG = {{1,3},{2,4}} and a non-zero prediction.**

### Row (f) numerical comparison — TENSOR case only

Predicted: **ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²)** at Reading B, V_TRUNC = 16.

Compute script: `C:/Collatz/d2_extract_phi_X_squared.py` (sandbox-denied this session; user must run).
Output target: `C:/Collatz/experiments_output/d2_phi_X_squared.json`.

**Per scalar reduction, comparing measured M_4_alt vs predicted product:**

| Reduction | Measured M_4_alt | Predicted ϕ(X1²)·ϕ(X2²) | Ratio | Match? |
|---|---|---|---|---|
| tr_pi | 5.357×10⁻² | _pending compute_ | — | — |
| vac_pi | 4.775×10⁻³ | _pending compute_ | — | — |
| delta_1 | 5.742×10⁻² | _pending compute_ | — | — |
| sum_entries | 6.089×10⁻¹ | _pending compute_ | — | — |

### What the row (f) numerical check decides

Three logical outcomes:

1. **EXACT MATCH** (ratio = 1.000 to ~10⁻³ at sum_entries): tensor-digraph BMT (or bigraph) fits row (f) quantitatively. But row (d) still mismatches — so even this is a PARTIAL fit only.
2. **QUALITATIVE MATCH, QUANTITATIVE MISMATCH** (predicted non-zero but ≠ measured at ratio ≠ 1): bigraph/BMT predicts the right STRUCTURE (non-zero crossing pair partition) but wrong magnitude. Confirms BMT moment formula is too loose for Syracuse's specific phase structure.
3. **PREDICTED VALUE NEAR ZERO** (ϕ(X̃_j²) ≈ 0 in Reading B): wouldn't expect this — second moments of centered operators are positive — but worth checking.

**Pre-numerical expectation:** outcome (2). Reading B's X̃_{j_2} has b_prior-dependent phase that couples non-trivially to the j_1 within-pair split. The measured 6.089×10⁻¹ at sum_entries equals 9.6 × ϕ(X̃_2 X̃_1 X̃_2)_average ≈ 9.6 × 6.347×10⁻² (Fubini inner from D1 fubini_grid_pairs). This is NOT obviously the product of two single-index second moments — those don't include the cross-step phase coupling.

---

## 5. Files

- This file: `C:/Collatz/D2_BMT_BIGRAPH_COMPARISON.md`
- Companion: `D2_BMT_BIGRAPH_VERBATIM.md`, `D2_BMT_BIGRAPH_PREDICTIONS.md`, `D2_BMT_BIGRAPH_DISPOSITION.md`
- Syracuse data: `C:/Collatz/experiments_output/{monotone_diagnostic_n3.json, n4_alternating_diagnostic.json}`
- Pending: `d2_phi_X_squared.json` (run `python d2_extract_phi_X_squared.py`)
