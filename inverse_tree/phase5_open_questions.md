# Phase 5 — Open-question analysis

Followup on the 5 questions from `inverse_tree_findings.md`. Source: `tree_d50.parquet`.

## Q1 — Within-class residue collapse

At d=50 (79,255 nodes), individual-level variance decomposition of branching (binary outcome):

- Grand-mean branching: 0.2636
- SS_total (Bernoulli N·p(1−p)): 15383.3
- SS_between (a★ classes): 242.0
- **η²(a★_6) = 0.0157**
- η²(a★_6 × r_6 cells, for reference): 0.0949

**Verdict:** 
a★_6 captures only 1.6% of branching variance — branching is NOT meaningfully aligned with a★_6. The Phase-2 'a★ stratifies branching' finding reflects pooled means, not individual-level structure.

**What is the right coarsening?** Branching is a deterministic function of n mod 6: branching ⟺ n ≡ 4 (mod 6) AND n > 4. So η² for a partition by n%6 should be ≈1.

- η² for a★_6 partition (7 classes): **0.0157**
- η² for n%6 partition (6 classes): **1.0000**
- η² for joint (a★_6, n%6) partition: **1.0000**

n%6 captures all branching variance trivially (it's the deterministic predicate). The Phase 2 "a★ stratifies branching" finding is real but indirect: each a★ class happens to mix even/odd residues in characteristically different proportions. The cleanest statement is that **a★_6 is correlated with parity composition, not with branching directly.**

Branching ratio by parity within each a★ class:

| a★_6 | n_even | br_even | n_odd | br_odd |
|---:|---:|---:|---:|---:|
| 1 | 19,452 | 0.3334 | 0 | n/a |
| 3 | 25,615 | 0.3323 | 5,144 | 0.0000 |
| 9 | 13,607 | 0.3344 | 6,735 | 0.0000 |
| 27 | 3,509 | 0.3374 | 3,634 | 0.0000 |
| 81 | 462 | 0.3030 | 933 | 0.0000 |
| 243 | 35 | 0.5143 | 121 | 0.0000 |
| 729 | 0 | n/a | 8 | 0.0000 |

Even-n rows show branching rates near 1/3 (the within-even rate when k % 3 cycles uniformly); odd-n rows show 0 (n odd ⟹ (n-1)/3 even ⟹ ineligible). a★_6 stratification of pooled branching is a re-expression of how each class's even/odd ratio varies.

Cell-level table: `phase5_q1_cells.parquet`. Parity table: `phase5_q1_parity.parquet`.

## Q2 — Source of the 0.264 branching asymptote

Decomposition over d=30..50:

- P(n ≡ 1 mod 3) in tree: **0.3333**  (uniform-integer baseline: 0.3333)
- P((n−1)/3 odd | n ≡ 1 mod 3) in tree: **0.7919**  (baseline: 0.5000)
- Product: **0.2639**
- Actual branching ratio: **0.2639**

Product matches actual to within 0.5%, so the 0.264 asymptote is fully explained by these two factors. Tree biases relative to integer-uniform: ΔP(mod 3 = 1) = -0.0000, ΔP(cond. eligible) = +0.2919. Empirical / naive (1/6) ratio: **1.584×**.

Per-depth decomposition: `phase5_q2_branching_decomp.parquet`.

## Q3 — Equilibrium from 7×7 transition matrix on a★_6

Built M[i,j] = E[# class-j children per class-i parent] from 300,345 parent–child edges (parents at d < 50).

- **λ_max = 1.2638** (predicted growth factor per layer)
- exp(empirical slope 0.234) = 1.2636

Predicted (left eigenvector) vs empirical equilibrium:

| a★_6 | empirical (d=50) | predicted (eigvec) | diff |
|---:|---:|---:|---:|
| 1 | 0.2454 | 0.2453 | -0.0002 |
| 3 | 0.3881 | 0.3884 | +0.0003 |
| 9 | 0.2567 | 0.2567 | +0.0001 |
| 27 | 0.0901 | 0.0899 | -0.0002 |
| 81 | 0.0176 | 0.0177 | +0.0001 |
| 243 | 0.0020 | 0.0019 | -0.0000 |
| 729 | 0.0001 | 0.0001 | -0.0000 |

Max |diff|: **0.0003**.

Eigenvector matches empirical equilibrium tightly — the 7-state Markov chain on a★_6 closes the loop.

Transition matrix: `phase5_q3_transition_matrix.csv`. Eigenvector table: `phase5_q3_eq_eigenvector.parquet`.

## Q4 — Rate of mod-2^k equidistribution

**Sample-size correction.** Raw chi² grows linearly with n_total even at fixed underlying distribution. The right metric is **chi²/n_total** (deviation per node).

Fit log(chi²/n_total) ~ a + b·d on window d=25..50:

| k | slope b | intercept a | half-life (layers) |
|---:|---:|---:|---:|
| 6 | +0.0222 | 0.587 | inf (no decay) |
| 8 | +0.0569 | -0.355 | inf (no decay) |
| 10 | +0.0957 | -1.722 | inf (no decay) |

Slope ~ 0 means residue distribution is essentially STATIC — neither converging to nor diverging from uniform within d≤50.
(Raw chi²/dof grows ~exp(0.234·d), tracking the layer-count growth rate, not residue convergence.)

Plot: `phase5_q4_chi2_decay.{svg,png}`. Per-(d,k) data: `phase5_q4_chi2_decay.parquet`.

## Q5 — Re-rooted self-similarity (multi-founder)

For each a★_6 class, take ALL founder nodes (n in class whose parent is not in same class), enumerate sub-tree, fit slope on last 20 offsets if sub-tree depth ≥ 10.

| a★_6 | n_founders | n_with_slope | slope_median | slope_MAD | max_sub_size |
|---:|---:|---:|---:|---:|---:|
| 1 | 19,452 | 1,879 | 0.2381 | 0.0403 | 186,434 |
| 3 | 50,211 | 4,837 | 0.2374 | 0.0387 | 379,596 |
| 9 | 45,937 | 4,433 | 0.2369 | 0.0391 | 379,598 |
| 27 | 20,641 | 1,985 | 0.2364 | 0.0396 | 379,600 |
| 81 | 4,953 | 451 | 0.2344 | 0.0424 | 11,168 |
| 243 | 637 | 49 | 0.2352 | 0.0434 | 892 |
| 729 | 34 | 2 | 0.1384 | 0.1384 | 80 |

Plot: `phase5_q5_rerooted_slopes.{svg,png}`. Per-class table: `phase5_q5_rerooted.parquet`.

**Verdict:** Median slopes within ~0.005 of the full-tree slope (0.2343) ratify self-similarity beyond the single-founder analysis. Per-class MAD shows tightness across founders within each class.