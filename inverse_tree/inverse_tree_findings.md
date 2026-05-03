# Inverse Collatz Tree — Findings (d ≤ 50)

Inverse tree rooted at n=1, built to depth 50. **379,600 nodes.** Each node annotated with prefix-decomposition fields (a★, c★, ℓ, α_det) at modular resolutions k ∈ {6, 8, 10}, computed via the existing `deterministic_prefix` utility in `experiments/01_alpha_decomposition.py`. OEIS A005186 cross-check passes through d=20; forward-σ verification passes on a 10,000-node random sample.

> **Note (post-Phase-5 revision):** Section 2 has been rewritten. The original Phase 2 reading "a★ stratifies branching" turned out to be a pooled-mean artifact — at the individual level, branching is determined exactly by n mod 6, and a★_6 captures only ~1.6% of branching variance. Original wording is preserved at the bottom in "Appendix: Phase 2 reading (superseded)". See Section 6 for the Phase 5 follow-up findings on all five originally-open questions.

## 1. Tree statistics summary

- Total nodes through d=50: **379,600**
- Layer-count growth slope (log L_d vs d, window d=10..50): **0.2343** (R² = 0.9998)
- Tao-style asymptote log(4/3) ≈ 0.2877 — finite-window slope is below this, consistent with pre-asymptotic behavior at d=50
- Branching ratio (fraction of nodes with both predecessors): **0.264** at d ≥ 20 (stable)
- Build wall-time: 0.072s; parquet on disk: 4.2 MiB

## 2. a★ vs branching — revised verdict

a★_6 partitions the 64 mod-64 residues into 7 buckets {1, 3, 9, 27, 81, 243, 729}.

**Density.** P(a★_6 | depth) equilibrates by d ≈ 25 to a stable distribution: a★=3 (38.8%), a★=9 (25.7%), a★=1 (24.5%), a★=27 (9.0%), a★=81 (1.8%), a★=243 (0.2%), a★=729 (0.01%). Three values carry 89% of the mass. This equilibrium is **explained exactly** by the 7×7 transition matrix on a★_6 (Section 6, Q3): leading left-eigenvector matches empirical to within 0.0003, and λ_max = 1.2638 = exp(0.234).

**Branching is NOT actually a★-stratified at the individual level.** Pooled-by-a★ branching ratios at d=50 do separate (0.333 → 0.100 across a★ ∈ {1, 3, 9, 27, 81}), but this is an artifact:

- η²(a★_6) for individual-level branching = **0.016** (1.6% of variance)
- η²(n mod 6) for individual-level branching = **1.000** (the deterministic predicate)

Branching is exactly equivalent to **n ≡ 4 (mod 6) AND n > 4**. Within every a★_6 class, even-n nodes branch at ~0.33 and odd-n never branch. The Phase 2 "stratification" was a re-expression of how each a★ class mixes even and odd r_6 values:

| a★_6 | n_even | br_even | n_odd | br_odd |
|---:|---:|---:|---:|---:|
| 1 | 19,452 | 0.3334 | 0 | n/a |
| 3 | 25,615 | 0.3323 | 5,144 | 0.0000 |
| 9 | 13,607 | 0.3344 | 6,735 | 0.0000 |
| 27 | 3,509 | 0.3374 | 3,634 | 0.0000 |
| 81 | 462 | 0.3030 | 933 | 0.0000 |

**Verdict:** a★_6 is correlated with parity composition, not with branching directly. The right coarsening for branching is mod-6 (or equivalently parity-of-r combined with n%3).

## 3. Self-similarity verdict

Sub-trees rooted at the smallest n in each a★_6 class:

| a★_6 | root_n | root_depth | sub_total | slope (last 20 offsets) | Δ vs full tree |
|---:|---:|---:|---:|---:|---:|
| 1 | 64 | 6 | 186,434 | 0.2341 | -0.0002 |
| 3 | 16 | 4 | 379,596 | 0.2340 | -0.0003 |
| 9 | 4 | 2 | 379,598 | 0.2340 | -0.0003 |
| 27 | 1 | 0 | 379,600 | 0.2340 | -0.0003 |
| 81 | 7 | 16 | 11,168 | 0.2328 | -0.0015 |

The "smallest n in class" methodology has a degeneracy: a★=27's smallest representative is n=1 itself, so its sub-tree is the whole tree; classes 1/3/9/27 cover ≥98.4% of the tree by construction. Phase 5 addresses this by enumerating ALL "founder" sub-trees per class (founders = nodes whose parent is not in the same class) and reporting median slope across many independent sub-trees:

| a★_6 | n_founders | n_with_slope | slope_median | slope_MAD |
|---:|---:|---:|---:|---:|
| 1 | 19,452 | 1,879 | 0.2381 | 0.040 |
| 3 | 50,211 | 4,837 | 0.2374 | 0.039 |
| 9 | 45,937 | 4,433 | 0.2369 | 0.039 |
| 27 | 20,641 | 1,985 | 0.2364 | 0.040 |
| 81 | 4,953 | 451 | 0.2344 | 0.042 |
| 243 | 637 | 49 | 0.2352 | 0.043 |

**Verdict:** Median slopes within ±0.005 of the full-tree slope (0.2343) across thousands of independent founder sub-trees. **Self-similarity is ratified.**

## 4. Surprises

1. **Branching asymptote 0.264 fully explained by conditional eligibility bias.** P(n ≡ 1 mod 3 | tree node) is unbiased at exactly 0.333. The bias is in the conditional: P((n−1)/3 odd | n ≡ 1 mod 3, in tree) = **0.792** (vs 0.500 baseline) — a 1.58× over-representation of n ≡ 4 (mod 6) over n ≡ 1 (mod 6). Product 0.333 × 0.792 = 0.264 matches the empirical rate exactly.

2. **a★_6 density converges fast and the closed form is exact.** By d=25 (227 nodes total at that layer), the distribution is within ~2% of its d=50 shape. Phase 5 derives the limit: the leading left-eigenvector of the 7×7 transition matrix on a★_6 matches empirical to **0.0003** across all 7 classes.

3. **Residue equidistribution does NOT happen.** Sample-size-corrected chi²/n_total has slope **+0.022 (k=6), +0.057 (k=8), +0.096 (k=10)** on log scale — slowly INCREASING with depth, not decreasing. Residue distributions are converging to non-uniform stationaries, just like a★_6 does.

4. **Pooled-by-a★ branching ratios LOOKED like a geometric sequence (0.83, 0.81, 0.74, 0.60 ratios) — but this is an artifact.** Each a★ class has a characteristic even/odd split. As a★ grows, the class's even-fraction drops (1.00 → 0.83 → 0.67 → 0.49 → 0.33 for a★ ∈ {1, 3, 9, 27, 81}), and pooled branching = even-fraction × ~0.33. The "geometric sequence" tracks the even-fraction sequence, not anything intrinsic to branching.

## 5. Phase 5 follow-up: answers to the originally-open questions

Each open question from the Phase 4 draft is now answered. Detail in `phase5_open_questions.md`.

**Q1 — Within-class residue collapse:** **Refuted.** η²(a★_6) = 0.016. a★_6 captures less than 2% of branching variance. The right partition is n%6.

**Q2 — Source of 0.264 branching:** **Solved.** P(n%3==1) is unbiased; P((n−1)/3 odd | n%3==1) = 0.792 vs 0.500 baseline. The full 1.58× excess over 1/6 lives in the conditional eligibility.

**Q3 — Closed-form a★_6 equilibrium:** **Solved.** 7×7 transition matrix M[i,j] = E[# class-j children per class-i parent] has λ_max = 1.2638 (= exp(0.234)) and a left-eigenvector matching empirical equilibrium to 0.0003. `phase5_q3_transition_matrix.csv` ships the matrix.

**Q4 — Mod-2^k equidistribution rate:** **Inverted.** Sample-size-corrected metric is essentially flat or slightly INCREASING — no convergence to uniform. Residues equilibrate to a non-uniform stationary.

**Q5 — Re-rooted self-similarity:** **Confirmed.** Multi-founder analysis across 6/7 classes: median slopes in [0.234, 0.238] across thousands of founder sub-trees per class. Within 0.005 of full tree.

## 5b. New open questions (raised by Phase 5)

1. **Why is P((n−1)/3 odd | n ≡ 1 mod 3) = 0.792 in the tree?** This is the entire source of the branching bias. The inverse-Collatz step that creates n ≡ 1 mod 3 nodes is "child = (n_parent − 1)/3", and that child is odd. So all n%3==1 nodes that arrived via the (n−1)/3 step are ODD. But many n%3==1 nodes also arrive via the doubling step. The 0.792 is the equilibrium odd-fraction among n%3==1 tree nodes. Predict from a 2-state chain on (n%3==1, parity) and verify.

2. **Closed form for the residue stationary distribution mod 2^k.** Q3 worked for the 7-state a★_6 chain. The same machinery should work at the full 64-state mod-64 level — derive the mod-64 transition matrix, take its leading eigenvector, compare to empirical residue distribution at d=50. If it matches, the slow chi²/n drift is explained as approach to that fixed point.

3. **Why does each a★_6 class have a characteristic even/odd split?** a★_6=1 is 100% even residues; a★_6=243 is 22.4% even. The deterministic_prefix function maps r → a_final independent of parity, so this asymmetry must reflect how the predecessor map distributes residues into a★ buckets. Worth a short structural derivation.

## 6. One-sentence-per-phase summary

**Phase 1:** Built and serialized the d≤50 inverse tree (379,600 nodes, 0.07s, OEIS PASS, forward-σ PASS, 4.2 MiB parquet). **Phase 2:** Computed branching, density, and self-similarity diagnostics; pooled-by-a★ branching ratios appeared stratified (later revised). **Phase 3:** Three figures emitted — tree to d=15 colored by a★_6, stacked-area density showing equilibrium by d≈25, per-stratum branching chart. **Phase 4:** Consolidated report (this file). **Phase 5:** All five originally-open questions resolved: a★_6 captures only 1.6% of branching variance (the apparent stratification was a parity-mix artifact); the 0.264 asymptote is fully explained by P((n−1)/3 odd | n%3==1) = 0.792; a 7×7 transition matrix on a★_6 reproduces the equilibrium distribution exactly (λ_max = 1.264, max eq diff = 0.0003); residue equidistribution is NOT happening — chi²/n is flat or slightly growing; multi-founder re-rooted self-similarity holds across thousands of independent sub-trees.

## 7. Files in `C:\Collatz\inverse_tree\`

- `tree_d50.parquet` — full annotated tree (4.2 MiB, 20 columns)
- `build_tree.py` / `phase2_diagnostics.py` / `phase3_visualize.py` / `phase5_open_questions.py` — runnable artifacts
- `phase1_checkpoint.md` / `phase2_findings.md` / `phase5_open_questions.md` / `inverse_tree_findings.md` — staged reports
- `phase2_*.parquet` — Phase-2 diagnostic tables
- `phase3_*.{svg,png}` — Phase-3 figures (tree, stacked-area density, branching-by-stratum)
- `phase5_q1_cells.parquet` / `phase5_q1_parity.parquet` / `phase5_q2_branching_decomp.parquet` / `phase5_q3_transition_matrix.csv` / `phase5_q3_eq_eigenvector.parquet` / `phase5_q4_chi2_decay.parquet` / `phase5_q5_rerooted.parquet` — Phase-5 diagnostic tables
- `phase5_q4_chi2_decay.{svg,png}` / `phase5_q5_rerooted_slopes.{svg,png}` — Phase-5 figures

## Appendix: Phase 2 reading (superseded by Phase 5)

The original Phase 2 verdict, retained for traceability:

> **Branching ratio is a★-stratified, not collapsed.** At d=50: 0.333 (a★=1), 0.277 (a★=3), 0.224 (a★=9), 0.166 (a★=27), 0.100 (a★=81), 0.115 (a★=243; small-sample). Mean spread across classes (d=30..50, cells with ≥10 nodes): 0.254. Mean within-depth std: 0.089.
>
> Verdict: a★_6 stratifies branching empirically — lines fan out rather than collapse.

The pooled-mean observation is correct as a description of the stratified means, but the η²=0.016 individual-level result shows the underlying mechanism is the per-class even/odd composition, not anything specifically about a★ encoding branching-relevant structure.
