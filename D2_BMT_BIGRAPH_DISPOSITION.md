# D2 — BMT / Bigraph Disposition (Fit / Partial / No fit)

**Date:** 2026-05-15
**Mode:** E.
**Companion:** `D2_BMT_BIGRAPH_VERBATIM.md` (verbatim defns), `D2_BMT_BIGRAPH_PREDICTIONS.md` (derivations), `D2_BMT_BIGRAPH_COMPARISON.md` (numerics).

---

## Verdict: **NO FIT** for BMT or bigraph independence on a 2-vertex Syracuse graph

Specifically:

| Row | Syracuse measured (sum_entries) | BMT (best digraph) | Bigraph (best pairwise) | Match? |
|---|---|---|---|---|
| (b) | ≈ 0 (noise floor 10⁻⁷) | 0 (every digraph) | 0 (every pairwise) | ✓ |
| (d) | 1.078×10⁻¹ | **0** (every digraph) | **0** (every pairwise) | ✗ DECISIVE |
| (f) | 6.089×10⁻¹ | 0 (mono) or ϕ(X̃_1²)·ϕ(X̃_2²) (tensor) | same as BMT | qualitatively under tensor only; quantitative pending |

**The row (d) mismatch is structural, exact, and impossible to escape** by any choice of digraph (BMT) or bigraph (Gilliers–Jekel). Both frameworks predict 0 for row (d) under every pairwise relation. Syracuse measures 0.108 (5–4 orders of magnitude above noise).

---

## 1. Sharp reasoning — why row (d) is decisive

**BMT Defn 3.4** (p9): ϕ(a_1 a_2 a_3) = ∏_{B ∈ ker_G[i]} ϕ[(a_k)|_B]. For i = (j_1, j_2, j_1) the position 2 with label j_2 is structurally a **singleton** in ker_G[i] regardless of digraph choice, because position 2 has no "twin" with label j_2 elsewhere in the word. So one block is {2} with ϕ((a_k)|_{2}) = ϕ(X̃_{j_2}) = 0 by Reading B centering. The product of all blocks is 0.

**Bigraph Defn 1.4** (p4): φ(a_1 ··· a_k) = ∑_{π ∈ P(c, G)} K^free_π(a_1, ..., a_k). For c = (j_1, j_2, j_1), every admissible partition π either keeps position 2 a singleton {2} (forcing κ_1^free(X̃_{j_2}) = 0 in the product), or fails Defn 1.2 condition (1) (same-block ⇒ same color, blocking partitions that bundle position 2 with positions of different label). Either way the sum vanishes.

Both formulas vanish at row (d) **structurally**, not just under the wrong digraph choice. There is no "loophole" digraph or bigraph that lifts row (d) off zero.

## 2. The auditor's qualitative claim at row (f)

The auditor's predictive check claimed: BMT monotone-digraph at row (f) predicts ker_G[(j_1,j_2,j_1,j_2)] = {{1,3},{2,4}} → ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²), non-zero, qualitatively matching Syracuse's 6.089×10⁻¹.

**Direct check from BMT Defn 2.8 + Prop 3.5(iii):** The "monotone digraph" has E = {(i, j) : j < i}. So **(j_2, j_1) ∈ E** (since j_1 < j_2) but **(j_1, j_2) ∉ E**.

- ker_G[i]: 1 ∼ 3 requires (i_2, i_1) = (j_2, j_1) ∈ E. **TRUE.** OK.
- 2 ∼ 4 requires (i_3, i_2) = (j_1, j_2) ∈ E. **FALSE.**
- So ker_G[i] = {{1, 3}, **{2}**, **{4}**}, NOT {{1,3},{2,4}}.

BMT prediction under monotone digraph = ϕ(X̃_1²) · ϕ(X̃_2) · ϕ(X̃_2) = ϕ(X̃_1²) · 0 · 0 = **0**, not non-zero.

**The auditor's claim was wrong.** Only the TENSOR (complete) digraph gives ker_G = {{1,3},{2,4}} and the non-zero ϕ(X̃_1²)·ϕ(X̃_2²) prediction.

## 3. What still hangs on the numerical check

Even with the auditor's specific BMT prediction corrected: BMT under the **tensor digraph** (and bigraph with E_2 ⊇ tensor edge between j_1, j_2) does predict ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²) ≠ 0 for row (f). Whether this magnitude matches 6.089×10⁻¹ requires the numerical computation of ϕ(X̃_1²) and ϕ(X̃_2²) under Reading B at V_TRUNC = 16.

Script `C:/Collatz/d2_extract_phi_X_squared.py` is ready; sandbox denied execution this session. **User action requested:** `python C:/Collatz/d2_extract_phi_X_squared.py` to populate `experiments_output/d2_phi_X_squared.json`.

But even an EXACT row (f) numerical fit under tensor BMT/bigraph would leave row (d) unfit. Tensor independence between j_1 and j_2 ALSO predicts row (d) = 0 (since the middle position 2 with label j_2 has no other position with label j_2 to pair with — singleton block under ker[i] regardless of digraph).

## 4. Sharper finding — Syracuse fits NO 2-vertex pairwise relation

The D1 numerical data already established Syracuse violates monotone (peak rule fails at row (d) for monotone, fails at row (f) for everyone). The new finding here is:

**Syracuse violates monotone, anti-monotone, Boolean, free, AND tensor on a 2-vertex graph between A_{j_1} = B⟨X̃_{j_1}⟩ and A_{j_2} = B⟨X̃_{j_2}⟩, in BOTH the BMT and bigraph frameworks.** Row (d) ≠ 0 is incompatible with every one of these.

The structural origin: row (d) ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) non-vanishing IS the Δ_{j_2}(b_{[1, j_1]}) phase twist from the **conditional-expectation level grading** (per `H1_PRIME_LOW_ORDER_CHECKS.md` §3.1). This phase twist is a B_marginal-measurable object encoding cross-step phase coupling that BMT and bigraph's formulas DON'T see because they evaluate ϕ block-wise on each kernel-block, and Reading-B ϕ(X̃_{j_2}) = 0 in isolation.

**Syracuse needs a framework where the conditional expectation ϕ at the peak position itself depends on the surrounding word context** — i.e. the level-graded reading ii of `H1_PRIME_LOW_ORDER_CHECKS.md` §3.1.a. BMT and bigraph are single-state moment relations and do not encode this.

---

## 5. Implication for the D2 audit

The D2 Tier 1 audit's "Claim 5 FAILED" verdict identified a real loophole (BMT/bigraph defns are moment relations not requiring iid copies), but **the loophole closes here**: even though BMT and bigraph DON'T require iid copies in the definition, the moment formulas they prescribe still don't fit Syracuse — they ALL predict row (d) = 0, contradicting measured 0.108.

Three rewrites of the audit verdict are appropriate:

- **D2_TIER1_DISPOSITION.md "Claim 5"** should be re-stated: BMT/bigraph DON'T require iid copies in the DEFINITION (correctly noted in the audit), but they STILL fail Syracuse at row (d) STRUCTURALLY (the kernel/admissible-partition condition isolates the peak position 2, forcing ϕ(X̃_{j_2}) = 0 to kill the prediction).

- The "Tier 2/3 loophole" the D2 audit flagged is **closed**: BMT and bigraph are not viable fits for Syracuse.

- The auditor's qualitative claim that BMT-monotone-digraph predicts ϕ(X̃_{j_1}²)·ϕ(X̃_{j_2}²) at row (f) is **wrong**; only BMT-tensor or bigraph-tensor give that prediction, and even then row (d) still mismatches.

---

## 6. What this means for downstream

- **c = 7/45 derivation:** unaffected. Per `D1_DISPOSITION.md` and `D3_DERIVATION_AUDIT.md`, the c=7/45 leading coefficient never depended on row (d) or row (f), so this disposition has zero impact on it.

- **Tier-N search for a fitting framework:** must continue. BMT and bigraph are eliminated as candidates. The level-graded conditional-expectation reading (H1_PRIME §3.1.a (ii)) remains the structural target; it requires a multi-state or filtered framework, which neither BMT nor bigraph provides.

- **Monotone-terminal verdict** (`OBSTRUCTION_MAP_TERMINAL.md`): Muraki 2003 monotone B-valued independence (with appropriate filtration) remains the load-bearing candidate per the prior MEMORY entry. BMT/bigraph being eliminated here REINFORCES that the right home is monotone-style with filtration, not a graph-mixture extension.

---

## 7. Files

- This disposition: `C:/Collatz/D2_BMT_BIGRAPH_DISPOSITION.md`
- Verbatim defns: `C:/Collatz/D2_BMT_BIGRAPH_VERBATIM.md`
- Predictions derivation: `C:/Collatz/D2_BMT_BIGRAPH_PREDICTIONS.md`
- Comparison table: `C:/Collatz/D2_BMT_BIGRAPH_COMPARISON.md`
- Numerical script (pending user run): `C:/Collatz/d2_extract_phi_X_squared.py` → `experiments_output/d2_phi_X_squared.json`

## 8. Mode-E gaps

- ϕ(X̃_{j_1}²), ϕ(X̃_{j_2}²) numerical values pending user run of `d2_extract_phi_X_squared.py` (sandbox denied this session). The verdict NO FIT is robust to this gap because row (d) decides; the row (f) numerical check would only refine the "qualitative match under tensor" line, not change the disposition.
- κ_2^free in the bigraph formula for cross-step variables j_1 ≠ j_2 (the (b) row in bigraph): assumed = ϕ(X̃_{j_1} X̃_{j_2}) by the standard NS06 partitioned free cumulant convention; this is exact for centered variables and gives the same value Syracuse measures at row (b) (≈ 0 by independence of within-pair splits). No quantitative gap here.
- The bigraph defn uses K^free_π for arbitrary π; for crossing pair partition {{1,3},{2,4}} and a single algebra A_{j_1} the NS06 partitioned free cumulant is the product κ_2(X̃_{j_1}, X̃_{j_1}) · κ_2(X̃_{j_2}, X̃_{j_2}) which for centered variables equals ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²). This identification is standard NS06 Chapter 11; explicit citation is in NS06 not pulled here.
