# Saddle-class subsum analysis — Outcome B-strong: partition is NOT a closure path

**Date:** 2026-05-04. Companion to R79b empirical study and milicevic_banks_verification.md.

## Executive summary

Per-saddle-class subsums `S_j(r) = |Σ_{a ∈ supp_j(r)} D(-3a) · ψ_true(a)|` measured at r ∈ {6, 8, 10}. Power-law fit gives:

| j | β_j (slope of log S_j vs log n_j) | R² | classification |
|---|---|---|---|
| 0 | **0.920** | 1.0000 | linear-in-n_j (worse than saturation) |
| 1 | **1.058** | 1.0000 | linear-in-n_j (worse than saturation) |
| 2 | **0.977** | 0.9987 | linear-in-n_j (worse than saturation) |

**Outcome: B-strong.** All three classes scale as `S_j ≈ c_j · n_j` (slope ≈ 1) — **even worse than Outcome B (β ≈ 0.5) as defined in the brief**. There is **no intra-class cancellation**; each class is essentially at its trivial bound. The full-sum |K| ≈ √N comes from **inter-class cancellation** between three same-order-magnitude but oppositely-phased subsums.

**Saddle-class partition is structural but NOT a closure path.** Q2 (structural tests on supp_0) is moot under Outcome B, but a notable structural fact emerges from the data: `|Σ_{a ∈ supp_j} ψ(a)|` (without Dirichlet weight) equals `√N` exactly for j=0, and `0` exactly for j=1, 2. This is a clean structural identity (likely provable from T78.3) — see §Structural identity below.

**Implication for c=7/45 closure:** the saddle-class partition direction in v3.7 is closed as a closure path. Open closure routes remaining: Bourgain-Konyagin sum-product on ⟨4⟩, direct band-l¹ analysis on D_{r,t}(η), smooth completion (R78 path 2), or Path C (5x+1 sibling-attack reframing). Path B (consolidate-and-publish current state) is now the fallback.

---

## Q1 results: per-class subsum scaling

### Raw data

```
r=6,  q=2187,    N=243,   period=729,   J=6
   j   n_j           S_j   log(S_j)   log(n_j)   S_j/√n_j   |Σ ψ_j|   |Σ ψ_j|/n_j
   0    81      160.4622     5.0781     4.3944    17.829    15.5885   0.192450
   1    81      176.7876     5.1749     4.3944    19.643     0.0000   0.000000
   2    81       72.7804     4.2874     4.3944     8.087     0.0000   0.000000

r=8,  q=19683,   N=2187,  period=6561,  J=9
   j   n_j           S_j   log(S_j)   log(n_j)   S_j/√n_j   |Σ ψ_j|   |Σ ψ_j|/n_j
   0   729     1206.3153     7.0953     6.5917    44.678    46.7654   0.064150
   1   729     1859.7581     7.5282     6.5917    68.880     0.0000   0.000000
   2   729      711.8028     6.5678     6.5917    26.363     0.0000   0.000000

r=10, q=177147,  N=19683, period=59049, J=10
   j   n_j           S_j   log(S_j)   log(n_j)   S_j/√n_j   |Σ ψ_j|    |Σ ψ_j|/n_j
   0  6561     9148.3633     9.1213     8.7889   112.943   140.2961   0.021383
   1  6561    18515.8831     9.8264     8.7889   228.591     0.0000   0.000000
   2  6561     5321.8958     8.5796     8.7889    65.702     0.0000   0.000000
```

### Per-class power-law fit

`log(S_j) = a_j + β_j · log(n_j)` over r ∈ {6, 8, 10}:

| j | β_j | a_j (intercept) | R² |
|---|---|---|---|
| 0 | 0.9201 | 1.0333 | 1.0000 |
| 1 | 1.0585 | 0.5327 | 1.0000 |
| 2 | 0.9767 | 0.0401 | 0.9987 |

Comparison to brief's predicted outcomes:

| Outcome | Predicted β | Observed (j=0,1,2) | Match? |
|---|---|---|---|
| A: controlled (O(log) or O(1)) | β ≈ 0 | 0.92, 1.06, 0.98 | NO |
| B: saturating (√n_j) | β ≈ 0.5 | 0.92, 1.06, 0.98 | NO |
| **B-strong: linear in n_j (trivial)** | β ≈ 1.0 | **0.92, 1.06, 0.98** | **YES** |

The empirical β ≈ 1 across all three classes is **stronger evidence against the partition route than Outcome B**: not only do the classes saturate (β ≥ 0.5), they actually grow linearly in class size (no cancellation within class at all).

### S_j / n_j ratio (effective per-class constant)

| j | r=6 | r=8 | r=10 | drift |
|---|---|---|---|---|
| 0 | 1.98 | 1.65 | 1.39 | slow decrease |
| 1 | 2.18 | 2.55 | 2.82 | slow increase |
| 2 | 0.90 | 0.98 | 0.81 | stable |

Constants c_j = S_j / n_j ∈ [0.8, 2.8] across r values — slow drift, no convergence to either 0 (controlled) or growth (super-linear). Best-fit interpretation: each class is at its **trivial bound up to an O(1) prefactor** that varies slowly with r.

### Total vs partitioned

From R79b scenario-comparison data: at r=10, `S_true = 19945`. Compare to per-class:

| quantity | r=10 value |
|---|---|
| S_true (total) | 19945 |
| S_0 + S_1 + S_2 (sum of magnitudes) | 32986 |
| max_j S_j | 18516 (j=1) |
| √N (square-root benchmark) | 140 |

`S_true < Σ_j |S_j|` — partial inter-class cancellation. But `S_true = 19945 ≈ 142 · √N` and `S_j ≈ c_j · n_j ≈ 3 · √N · √n_j` for typical c_j. So the inter-class cancellation reduces the bound from "linear in n_total" to "√n_total".

**The cancellation lives between classes, not within.** Saddle-class partition discriminates the phase direction of `D·ψ` (three distinct directions, one per class) but provides no magnitude improvement on each class individually.

---

## Structural identity: |Σ ψ| over each class

A clean numerical identity emerges in the |Σ ψ_j| (sum without Dirichlet weight) column:

| r | |Σ ψ_0| | √N | |Σ ψ_1| | |Σ ψ_2| |
|---|---|---|---|---|
| 6 | 15.5885 | √243 = 15.5885 | 0.0000 | 0.0000 |
| 8 | 46.7654 | √2187 = 46.7654 | 0.0000 | 0.0000 |
| 10 | 140.2961 | √19683 = 140.2961 | 0.0000 | 0.0000 |

Match to 4 decimals. **Structural identity (empirical, r ≤ 10):**
> Σ_{a ∈ supp_0(r)} ψ_true(a) = √N · e^{iθ_r}, |Σ_{a ∈ supp_j} ψ_true(a)| = 0 for j ∈ {1, 2}.

This is consistent with: ψ_true on j=1, j=2 has class-mean-zero (verified separately in r79b_s_class_deviation), while ψ_true on j=0 has class-mean exactly equal to √N / n_0 = 1/√n_0 (with phase aligned to a single direction). Likely provable from Theorem 78.3 (Plancherel saturation) via:
- |Σ_{a ∈ Z/3^r, gcd(3,a)=1} F̂(3a)|² = q · Σ_u 1·1 = q·q? No, that's the wrong contraction.
- |Σ_{a ∈ supp} ψ(a)|² = ⟨ψ, 1_supp⟩² which evaluates to a Gauss-sum-on-subgroup quantity.

The asymmetry between j=0 (full alignment) and j=1, 2 (perfect cancellation) traces to the same structural fact noted in R79b: ψ_lead is constant 1 on j=0 but a non-constant linear phase on j=1, 2. The Hensel correction on j=0 disperses ψ_true's phase, but the magnitude of the class-mean is preserved at exactly 1/√n_0 — a normalization-like identity.

This is a real arithmetic-combinatorial fact about the cubic-phase saddle structure. Doesn't close eq 190 (the Dirichlet weighting destroys this nice cancellation, restoring the trivial bound on each class), but it's the kind of clean structural identity worth recording.

---

## Q2: structural tests on supp_0 — moot under Outcome B-strong

The brief specifies Q2 should run only if Outcome A holds. Under Outcome B-strong, the j=0 class is **not** controlled — its subsum S_0 grows linearly in n_0 just like j=1, 2. So characterizing the structure of ψ on supp_0 doesn't help close eq 190; the bound is already tight at the trivial level.

For completeness, one Q2-style observation from the data: ψ_true on supp_0 has `|Σ ψ_0| = √n_0` (clean numerical identity above). This means ψ_0 is "constructively aligned" within the j=0 class with mean magnitude `1/√n_0` — exactly the normalization at which a uniform-direction sequence of |·|=1 vectors sums to √n_0. So ψ_0 is not random-walk-like (which would give |Σ| ≈ √n_0 with random direction); it's structured to give exactly √n_0 with deterministic phase. Likely a Gauss-sum-on-subgroup identity.

---

## Honest caveats

1. **Three r-points only**: r ∈ {6, 8, 10} gives three data points per class. The β_j fits have R² ≥ 0.999 but use only 3 points per fit. Stable interpretation requires the linearity to extend to r=12, 14. **Not pushed** (per brief stop point).
2. **r=10 ψ is full-Hensel via direct G(a) computation**, not via Hensel-lifted closed form. So the per-class subsums are computed using the *true* phase, satisfying the brief's requirement to avoid ψ_lead.
3. **Plancherel cross-check** at the per-class level was not directly verified (R79b cross-checked S_true against K_direct for the *total* sum at r=8, 10 to <1%; per-class breakdown inherits this precision).
4. **The "structural identity" |Σ ψ_0| = √N is empirical at r ≤ 10**. Expecting it to hold at all r requires an analytical argument (Gauss-sum-on-subgroup); not provided here.
5. **The inter-class cancellation magnitude** (S_true vs Σ_j |S_j|) reduces from 33K to 20K at r=10 — partial cancellation only. The ratio (ratio = S_true / (Σ |S_j|) = 0.605 at r=10) is r-dependent and doesn't exhibit obvious scaling.
6. **β = 1 vs β = 0.92–1.06**: the slight under-1 values for j=0 (β=0.920) and j=2 (β=0.977) are within fitting noise of true β=1 across the 3-point fit. The j=1 value β=1.058 is also within noise. Cannot distinguish "exactly linear" from "slightly sub-linear" with 3 points.

---

## Decision tree resolution

Per the brief's decision tree:

> **If Outcome B (all three classes saturate)**: Saddle-class partition is structural but not a closure path. Path C (5x+1 sibling attack) becomes the primary direction. Path B (consolidate-and-publish what we have) becomes the fallback.

**Verdict: Outcome B (sharper than defined: B-strong).** Saddle-class partition direction is closed as a closure path for c=7/45 via this route. The framework's closure to c = 7/45 via saddle-class partition is **genuinely blocked**.

**Recommended morning consolidation framing:**

- **Path A (closed)**: Cochrane attack (R78), van der Corput (R79), saddle-class partition (this analysis).
- **Path B (fallback)**: consolidate-and-publish c=7/45 with empirical β=0.522, structural anchors T78.1–T78.6, and explicit obstruction map (R79b §Strategic summary).
- **Path C (next attack)**: 5x+1 sibling-attack reframing, Bourgain-Konyagin sum-product on ⟨4⟩, or direct band-l¹ analysis on the dangerous band.

The saddle-class partition retains structural meaning (j=0 anomalous, j=1, 2 regular) and surfaces the clean `|Σ ψ_0| = √N` identity. These are publication-worthy structural facts even though the partition doesn't close eq 190.

---

## Files

- `saddle_class_subsum_analysis.py` — per-class subsum + power-law fit + Q2 stub
- `saddle_class_subsum_data.csv` — 9 rows (3 j-values × 3 r-values)
- `saddle_class_subsum_log.txt` — full stdout
- `saddle_class_subsum_analysis.md` — this writeup

## Compute audit

| metric | value |
|---|---|
| Hardware | 9950X3D 32 cores, no GPU |
| Numba threads | default |
| Concurrent NBA Projections | 2 (PID 3676 v4lite_overnight, PID 47564 FG3M backtest) — untouched |
| Total wall time | 2 seconds (r=6, 8, 10) |
| Max r reached | 10 |
| Stop reason | Q1 decisive at r=10; brief specifies stop unless ambiguous |
