# Phase 2 Diagnostics — Inverse Collatz Tree (d ≤ 50)

Source: `tree_d50.parquet` (379,600 nodes).

## 2.1 Branching ratio by depth

| depth | n_nodes | n_branching | branching_ratio |
|---:|---:|---:|---:|
| 5 | 2 | 0 | 0.0000 |
| 10 | 6 | 2 | 0.3333 |
| 15 | 24 | 5 | 0.2083 |
| 20 | 72 | 19 | 0.2639 |
| 25 | 227 | 60 | 0.2643 |
| 30 | 732 | 194 | 0.2650 |
| 35 | 2,365 | 623 | 0.2634 |
| 40 | 7,628 | 2,007 | 0.2631 |
| 45 | 24,561 | 6,464 | 0.2632 |
| 50 | 79,255 | 20,889 | 0.2636 |

Mean branching ratio over d=30..50: **0.2639**.

(Heuristic: ~1/3 of nodes are ≡1 mod 3; of those, roughly 1/2 yield odd (n−1)/3, giving an upper-bound expectation ~1/6 ≈ 0.167. Observed equilibrates around the empirical value above.)

## 2.2 a★_6 stratified branching — collapse test

Across d=30..50 (classes with ≥10 nodes): mean spread (max − min branching_ratio) = **0.2544**, mean std = **0.0886**.

**Verdict:** 
Branching ratio is a★-STRATIFIED — clear separation across classes indicates a★ encodes branching-relevant structure.

Per-class branching ratios at d=50 (≥10 nodes):

| a★_6 | n_nodes | branching_ratio |
|---:|---:|---:|
| 1 | 19,452 | 0.3334 |
| 3 | 30,759 | 0.2767 |
| 9 | 20,342 | 0.2237 |
| 27 | 7,143 | 0.1658 |
| 81 | 1,395 | 0.1004 |
| 243 | 156 | 0.1154 |

## 2.3 a★_6 density at deepest layer

P(a★_6 | depth=50):

| a★_6 | count | P |
|---:|---:|---:|
| 1 | 19,452 | 0.2454 |
| 3 | 30,759 | 0.3881 |
| 9 | 20,342 | 0.2567 |
| 27 | 7,143 | 0.0901 |
| 81 | 1,395 | 0.0176 |
| 243 | 156 | 0.0020 |
| 729 | 8 | 0.0001 |

Full table (all d): `phase2_astar_density.parquet`.

## 2.4 Self-similarity per a★_6 class

Sub-tree rooted at smallest n in each a★_6 class. Slope fit on last 20 offset levels. Full-tree slope (d=10..50) = **0.2343** for comparison.

| a★_6 | root_n | root_depth | sub_total | slope (last 20) | Δ vs full | branching_ratio_sub |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 64 | 6 | 186,434 | 0.2341 | -0.0002 | 0.2638 |
| 3 | 16 | 4 | 379,596 | 0.2340 | -0.0003 | 0.2638 |
| 9 | 4 | 2 | 379,598 | 0.2340 | -0.0003 | 0.2638 |
| 27 | 1 | 0 | 379,600 | 0.2340 | -0.0003 | 0.2638 |
| 81 | 7 | 16 | 11,168 | 0.2328 | -0.0015 | 0.2639 |
| 243 | 39 | 34 | 17 | 0.0000 | -0.2343 | 0.0000 |
| 729 | 127 | 46 | 7 | nan | +nan | 0.2857 |

**Verdict** (restricted to well-sampled sub-trees, ≥100 nodes — covers 5/7 classes, **100.00%** of total mass):

max |slope deviation vs full tree| = **0.0015**.
Sub-trees match the full tree's growth rate within tight tolerance — the tree is approximately self-similar with respect to a★_6 stratification across well-sampled classes.

Small-sample classes excluded from the verdict: 2 (0.00% of mass) — their sub-trees are too shallow inside d≤50 to fit a meaningful slope.

## 2.5 Residue-class density at d=50 — equidistribution

Chi² distance from uniform across observed residues (residues that don't appear are excluded from dof).

| k | M=2^k | unique residues observed | n_total | chi² | chi²/dof |
|---:|---:|---:|---:|---:|---:|
| 6 | 64 | 64 / 64 | 79,255 | 378088.2 | 6001.40 |
| 8 | 256 | 254 / 256 | 79,255 | 735945.3 | 2908.87 |
| 10 | 1024 | 909 / 1024 | 79,255 | 1229087.4 | 1353.62 |

Per-k full density tables: `phase2_residue_density_k{6,8,10}.parquet`.

Note: chi²/dof ≈ 1 is consistent with uniform; >> 1 indicates residue-dependent over- or under-representation.

## Outputs

- `phase2_branching.parquet` — depth, n_nodes, n_branching, branching_ratio
- `phase2_astar_branching.parquet` — (depth × a_final_6) cross-tab
- `phase2_astar_density.parquet` — long-format P(a★_6 | depth)
- `phase2_residue_density_k{6,8,10}.parquet` — long-format P(r_k | depth)
- `phase2_selfsimilarity.parquet` — sub-tree growth-rate per class