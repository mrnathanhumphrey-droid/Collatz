# Pre-Registration: F̂_p Theorem Candidate — Adversarial Verification

**Pre-registered:** 2026-05-11, before any compute on this task.
**Author:** Claude.
**Parent work:** Move 2 attempt (commit `45af179`, pre-reg `f96fb86`). Move 2 produced a candidate family-level theorem at the F̂_p level. This document locks the rules for adversarial verification of that candidate before any integration into Paper 4 or formalization as a standalone result.

---

## §0. The candidate theorem under verification (exactly as stated)

**THEOREM CANDIDATE (qx+1 Plancherel saturation, F̂ level):** For every prime `p ≥ 3` and every `r ≥ 2`, define

> `f_p(u) := e_M(c · (1+p)^u)`,  M = p^{r+1},  c ∈ (Z/M)^× (default c=1).

f_p has period p^r in u (since (1+p) has multiplicative order p^r in (Z/M)^× for any odd prime p ≥ 3 and r ≥ 1).

The **full-period Fourier transform** F̂_p^full of f_p (i.e., the DFT of f_p extended periodically from length p^r to length M) is supported on the principal-unit sub-support

> `supp(F̂_p^full) = { p·a (mod M) : a ∈ Z/p^r, a ≡ 1 (mod p) }`, with `|supp| = p^{r-1}`,

and on this support has uniform magnitude

> **|F̂_p^full(ξ)| = p^{(r+3)/2}**.

**Two normalizations (algebraically equivalent, both verified):**

| Object | Definition | Predicted magnitude on supp |
|---|---|---|
| **F̂_p^short(ξ)** | `Σ_{u=0}^{p^r − 1} f_p(u) · e_M(−ξu)` (length-period sum) | `p^{(r+1)/2}` |
| **F̂_p^full(ξ)** | `Σ_{u=0}^{M − 1} f_p^periodic(u) · e_M(−ξu)` (length-M extended) | `p^{(r+3)/2}` |

On the support `ξ = p·a`, the kernel `e_M(−p·a·u) = e_{p^r}(−au)` collapses, and `F̂_p^full(p·a) = p · F̂_p^short(p·a)` exactly (M/p^r = p copies of the period sum). The candidate theorem can be stated as `|F̂_short| = p^{(r+1)/2}` or `|F̂_full| = p^{(r+3)/2}`; **this verification tests both forms simultaneously, and the disposition is identical**.

**Empirical state (Move 2):** `qx1_move2_phase2_check.csv` verified the F̂_short magnitudes match p^{(r+1)/2} to ~1e-14 at six cells (p ∈ {3, 5, 7}, r ∈ {2, 3}). The CSV's `match=False` flag reflects a labeling mismatch (predicted column used the F̂_full normalization, computed values were F̂_short). This is documented in Move 2 attempt §A4. **The 6-cell empirical verification is real; this work extends and adversarially probes it.**

---

## §1. Hypotheses (locked)

**Primary hypothesis H_THEOREM:** The candidate magnitude formula

> `|F̂_p^short(ξ)| = p^{(r+1)/2}`  (equivalently `|F̂_p^full(ξ)| = p^{(r+3)/2}`)

holds **exactly** (to numerical precision) for every prime p ≥ 3 and every r ≥ 2, on the support `{p·a : a ≡ 1 mod p}` of cardinality p^{r-1}, with F̂ vanishing off this support. The proof template — Cochrane Theorem 2 + Plancherel + principal-unit Gauss-sum equidistribution — applies uniformly with no hidden p-dependence in any constant.

**Null hypothesis H_NULL:** At least one of the following:
- (N1) The candidate fails at some (p, r) outside the verified range — magnitude not equal to p^{(r+1)/2}, or not uniform on the predicted support, or non-zero values off the predicted support.
- (N2) The proof template has a p-dependent constant the family-level statement doesn't capture.
- (N3) The "exact equality" claim has finite-precision artifacts that break at higher r (the formula holds only approximately, not exactly).
- (N4) Support characterization (the `{a ≡ 1 mod p}` subgroup) needs refinement at some p.

**Pre-registered favoring NULL.** The candidate has been verified at 6 cells. Extension is the test. Conservative null-favored decision rules below.

---

## §2. Locked procedure (5 phases)

### Phase 1 — Extension to higher primes

**Test cells:** `p ∈ {11, 13, 17, 19, 23, 29, 31}` × `r ∈ {2, 3}` = 14 cells.

**Method:** Compute F̂_p^short at each cell using zero-padded FFT of length M (gives all ξ ∈ Z/M simultaneously). For each cell, report:
- Max magnitude on predicted support `{p·a : a ≡ 1 mod p}`.
- Min magnitude on predicted support.
- Max magnitude off predicted support (verify zero).
- Max relative deviation `|magnitude_observed − p^{(r+1)/2}| / p^{(r+1)/2}` on support.
- Support cardinality (verify = p^{r-1}).

**Decision rule:** Phase 1 PASSES if max relative deviation < 1e-12 (machine precision) across all 14 cells, AND off-support magnitudes < 1e-10 (effective zero), AND support cardinality matches predicted exactly. Phase 1 FAILS if deviation > 1e-8 at any cell, OR off-support magnitudes exceed 1e-8.

### Phase 2 — Extension to higher r

**Test cells:** `p ∈ {3, 5, 7}` × `r ∈ {4, 5, 6}` = 9 cells.

**Method:** Same as Phase 1. Higher r tests whether the formula's exactness holds at scale where periods reach p^6 ∈ {729, 15625, 117649}.

**Decision rule:** Same precision thresholds as Phase 1. Pass / fail per cell + per phase verdict.

**Memory note:** at p=7, r=6 the array length M = p^{r+1} = 7^7 = 823,543 — large but FFT-able. If memory becomes limiting, document the ceiling and stop at the highest feasible r per prime.

### Phase 3 — Proof template adversarial walkthrough at p = 11

The "p-blind proof template (Cochrane T2 + Plancherel + principal-unit equidistribution)" claim needs scrutiny at a prime not in the original verification set. **Walk through the proof at p = 11 explicitly** and verify:

1. **Cochrane Theorem 2** applies to the polynomial `g(u) = c·(1+11)^u mod 11^{r+1}`. Specifically: does the p-adic degree `D = deg_p H+` computation work the same way at p = 11 as at p = 3? Is the bound formula's dependence on p captured by the family-level statement?

2. **Plancherel identity** `Σ_ξ |F̂_p|² = M·Σ_u |f_p|²` holds at p = 11 with no modification (it holds at every finite cyclic group; no p-dependence). Confirm and document.

3. **Principal-unit equidistribution.** At p = 11, the principal units `{1 + 11k : k ∈ Z/11^r} ⊆ (Z/11^{r+1})^×` form a cyclic group of order 11^r (standard result, no p-dependence for odd primes). The Gauss-sum equidistribution argument over this cyclic group: at p = 11, does the same character-theoretic computation yield uniform magnitude on the sub-support?

For each of the three steps, the deliverable is a paragraph stating:
- The argument template at general p
- The specialization to p = 11
- Whether any constant in the argument depends on p in a way not captured by the family-level statement

**Decision rule:** Phase 3 PASSES if all three proof steps apply at p = 11 with no hidden p-dependence beyond what the family-level statement (uniform magnitude p^{(r+1)/2} or p^{(r+3)/2}) already absorbs. Phase 3 FAILS or PARTIALLY PASSES if any step has p-dependent constants not in the family-level claim.

### Phase 4 — Boundary behavior

**Boundary 1: p = 2.** The theorem excludes p = 2. Test why: compute F̂_2^short for (p=2, r ∈ {2, 3, 4}). Determine:
- Does (1+2) = 3 have full order 2^r in (Z/2^{r+1})^×?
- Does the predicted magnitude formula p^{(r+1)/2} hold at p = 2?
- If not, document the mode of failure.

**Boundary 2: r = 1.** The theorem requires r ≥ 2. Test r = 1 at p ∈ {3, 5, 7}:
- Compute F̂_p^short and check whether the predicted magnitude holds.
- Document why r = 1 is excluded (likely: support degenerates, equidistribution argument needs r ≥ 2).

**Boundary 3: Evaluation ceiling.** Identify the practical (p, r) ceiling beyond which the FFT array `(p^{r+1},)` becomes infeasible on current memory (~32 GB RAM). Document the ceiling.

**Decision rule:** Phase 4 produces documentation. If unexpected behavior at boundaries informs the theorem statement, surface in disposition.

### Phase 5 — Support characterization

The theorem statement says support = `{p·a : a ≡ 1 mod p}` of size p^{r-1}, with F̂ exactly zero outside. Verify by computing FULL F̂_p^short array for several cells and checking:

1. **Support cardinality:** `|{ξ ∈ Z/M : |F̂_p^short(ξ)| > 1e-10}|` matches predicted p^{r-1}.
2. **Support structure:** the numerical support is exactly the set `{p·a (mod M) : a ∈ Z/p^r, a ≡ 1 (mod p)}`. Compute the symmetric difference; should be empty.
3. **Off-support magnitudes:** max |F̂_p^short(ξ)| over ξ ∉ predicted support. Should be effective zero (< 1e-10).

**Test cells:** the 14 Phase-1 cells (cheap once FFT is computed). At higher-r cells from Phase 2, spot-check at least one cell per prime.

**Decision rule:** Phase 5 PASSES if support matches predicted set exactly across all tested cells. FAILS if predicted support is missing any element or contains extras. Symmetric-difference size > 0 indicates either implementation error or theorem-statement refinement needed.

---

## §3. Adversarial safeguards (locked)

### A1 — Numerical precision check (dual precision)

At three representative cells — (3, 2), (11, 2), (5, 3) — compute F̂_p^short via two independent methods:
- **Method 1:** numpy.fft.fft (float64).
- **Method 2:** mpmath direct summation at 50-digit precision.

Compare magnitudes on support and verify agreement to ≥ 1e-14. Disagreement beyond this indicates numerical artifacts; report explicitly. Agreement confirms the FFT-based bulk computation is faithful at the precision claimed.

### A2 — Implementation cross-check

Hand-compute |F̂_3^short(ξ=3)| at (p=3, r=2). Expected value: `√27 ≈ 5.196152422706632`. Hand-computation uses:
- f(u) = e^{2πi · k_u / 27} where k_u ∈ {1, 4, 16, 10, 13, 25, 19, 22, 7} are 4^u mod 27.
- F̂(3) = Σ e^{2πi(k_u − 3u)/27}, phases (k_u − 3u) mod 27.
- Six phases at 1/27, three at 10/27 → |F̂(3)|² = 27.

Verify the FFT implementation matches to 14 digits. (This was already verified implicitly in Move 2; A2 makes it explicit and traceable.)

### A3 — Honest deviation logging

Any phase that fails or produces unexpected output is documented precisely. No "rounding error" claims without precision-comparison evidence (A1). If a phase reveals scope-narrowing, the narrowing is named precisely (which p, which r, what deviation magnitude).

### A4 — Pre-reg adherence

Any deviation from this pre-registration during compute is logged in the results document with explicit reason. If a phase reveals the theorem needs reformulation, the reformulation is documented; the original claim is not silently revised.

---

## §4. Decision rules (locked, null-favored)

- **THEOREM_VERIFIED.** Phases 1–5 all pass. A1, A2 clean. Theorem candidate ranked as **verified up to the evaluation ceiling identified in Phase 4**, ready for formalization. The exact equality `|F̂_p^short(ξ)| = p^{(r+1)/2}` (equivalently `|F̂_p^full(ξ)| = p^{(r+3)/2}`) holds on predicted support, with predicted support structure, for all tested cells.

- **THEOREM_VERIFIED_WITH_SCOPE_NARROWING.** Some phases pass cleanly; others reveal the theorem holds with a refined scope. Examples:
  - Phase 3 reveals a p-dependent constant in proof template not captured by the family statement → theorem holds in the verified empirical range but proof template needs refinement.
  - Phase 4 reveals the practical ceiling sits earlier than expected → verified range is narrowed.
  - Phase 5 reveals support characterization slightly different at some primes.
  Document refined statement precisely.

- **THEOREM_FALSIFIED.** At least one Phase 1, 2, 4 (within boundary scope), or 5 cell produces clear failure: magnitude not equal to p^{(r+1)/2} on support, or non-zero off support beyond precision artifact, or support structure mismatched. Candidate retracted or reformulated; stop and report immediately.

- **INCONCLUSIVE.** Numerical precision issues prevent clean adjudication at higher cells (Phase 2 ceiling lower than r = 6 for all primes, A1 disagreement, etc.). Verified range stays at Phase 1's cells. Practical ceiling documented honestly.

If THEOREM_FALSIFIED triggers at any phase: stop subsequent phases. Don't waste compute on a falsified candidate.

If Phase 3 fails alone (other phases pass): disposition is THEOREM_VERIFIED_WITH_SCOPE_NARROWING (empirical pass, proof template needs refinement) or similar — the data still verifies but the rigor claim does not.

---

## §5. Deliverables

- `FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md` — this document, committed before compute.
- `FHAT_THEOREM_VERIFICATION_RESULTS.md` — disposition at top, all five phases documented, A1–A4 records, references to artifacts.
- `fhat_verification.py` — main verification script (FFT-based bulk computation, support / magnitude checks for Phases 1, 2, 4, 5).
- `fhat_verification_dual_precision.py` — A1 mpmath cross-check.
- `fhat_verification_results.csv` — per-cell magnitudes and deviation statistics for Phases 1, 2, 4, 5.
- `fhat_verification_a1_dual_precision.csv` — A1 mpmath comparison data.

---

## §6. Honest scope-of-attempt statement

This is verification of a candidate theorem produced in Move 2. The 6-cell empirical verification is real (CSV evidence). The candidate hasn't been formalized as a proven theorem yet; the proof template `Cochrane T2 + Plancherel + principal-unit equidistribution` is sketched in Move 2 attempt §Phase 2 Route (a) but not rigorously written.

This verification's purpose: **test the candidate adversarially before formalization**. If the candidate passes extension to 7 more primes (p up to 31) and 3 more r values (up to r = 6), and the proof template survives walkthrough at p = 11, and support / boundary / precision checks all pass, then the candidate is ranked as verified up to the evaluation ceiling. This is the correct epistemic standard before integrating into Paper 4 or proposing as a standalone result.

Conservative null-favoring is the discipline. If verification falsifies, that catch before publication is the discipline working — document and move on. If verification narrows scope, that refinement is itself a useful result.

Pre-registration locked.
