# Esscher-tilt closure for R58 residuals — outcome (γ); proposed closure REJECTED, residuals come from different mechanism

**Status.** Decisive. The Esscher-tilt-by-σ_orbit closure proposed for R58's
QSD-extreme residuals (Audit Claim 6 "propose-not-demonstrate") **does
not work**. Sweeping λ ∈ [-2, 2] gives best λ ≈ −0.01 (essentially zero),
with Pearson improvement only +0.014 (0.857 → 0.870). Train-test split
confirms identical result (no overfitting either way).

R58's residuals at r=5, r=23, r=13 are **NOT explained by σ-band /
depth-tilt structure**. They reflect a different mechanism — most likely
the structural distinction between inverse-tree subtree-size measure and
forward-trajectory measure (R69 / K vs P pattern).

**For audit hygiene:** convert Claim 6 from "propose-not-demonstrate" to
"REJECTED — Esscher tilt at σ_orbit doesn't close residuals." Walk back
the proposed closure honestly.

## Step 1 — Inverse tree at N=2^22 reproduced

| metric | value |
|---|---:|
| Tree size | 1,247,706 odd nodes |
| Max depth | 152 |
| Build time | 0.8s |
| **Baseline R58 Pearson** | **+0.8568** (vs paper 0.857 ✓) |
| Baseline MAE | 0.1189 |

Reproduces R58 paper's Pearson 0.857 to 4 decimals.

## Step 2 — Per-residue residuals confirm R58's QSD extremes

| r | D_R58 | D_emp | residual (D_emp − D_R58) |
|---:|---:|---:|---:|
| 1 | 1.240 | 1.217 | −0.022 |
| 3 | 0.972 | 0.990 | +0.018 |
| **5** | **1.376** | **1.836** | **+0.460** ← max enhancement under-predicted |
| 7 | 1.014 | 0.912 | −0.101 |
| 9 | 0.915 | 0.794 | −0.121 |
| 11 | 0.979 | 0.849 | −0.129 |
| **13** | **0.889** | **0.680** | **−0.209** ← max depletion under-predicted |
| 15 | 0.989 | 1.093 | +0.104 |
| 17 | 1.151 | 1.010 | −0.140 |
| 19 | 0.816 | 0.771 | −0.045 |
| 21 | 0.864 | 0.862 | −0.002 |
| **23** | **1.043** | **1.351** | **+0.308** ← enhancement under-predicted |
| 25 | 0.842 | 0.702 | −0.140 |
| 27 | 1.005 | 0.967 | −0.038 |
| 29 | 1.124 | 1.112 | −0.011 |
| 31 | 0.783 | 0.835 | +0.052 |

R58 paper's flagged extremes (r=5, r=23, r=13) confirmed as top 3 absolute
residuals.

## Step 3 — Esscher-tilt λ sweep: REJECTED

Tilted weight: w_λ(m) = subtree_size(m) · exp(λ · depth(m)).

| λ | Pearson | MAE |
|---:|---:|---:|
| −2.0 | very negative | high |
| −0.5 | +0.553 | 1.04 |
| **0.0 (baseline)** | **+0.857** | **0.119** |
| **−0.01 (optimal)** | **+0.870** | 0.118 |
| +0.5 | −0.209 | 0.683 |
| +1.0 | −0.256 | 1.08 |
| +2.0 | very negative | high |

**Best λ from refined sweep: −0.0100** (essentially zero).
**Best Pearson: +0.8704** (improvement only +0.014).

The Pearson surface has a sharp peak at λ ≈ 0 with rapid degradation in
either direction. Any meaningful tilt destroys the prediction.

## Step 4 — Per-residue at λ_optimal = −0.01

| r | baseline residual | tilted residual | change |
|---:|---:|---:|---:|
| 5 | +0.460 | +0.335 | improved by 0.13 |
| 13 | −0.209 | −0.225 | **worse by 0.02** |
| 23 | +0.308 | +0.288 | improved by 0.02 |
| 17 | −0.140 | −0.180 | worse by 0.04 |

The negative tilt slightly improves r=5 but slightly worsens r=13 and r=17.
**The proposed mechanism (Esscher tilt resolves QSD-tilt residuals) is
qualitatively wrong** — it can't simultaneously fix both enhancement
extremes (r=5, r=23) and depletion extreme (r=13).

## Step 5 — Train-test split confirms (no overfitting either way)

| metric | value |
|---|---:|
| λ_train (from tree at 2^21) | −0.0100 |
| Pearson_train | +0.8673 |
| Pearson_test (apply λ_train at 2^22) | +0.8704 |
| Improvement vs baseline | +0.0136 |

Train and test give the same answer. Not an overfitting issue. The Esscher
tilt simply doesn't help structurally.

## Step 6 — σ-quartile depth structure

| metric | value |
|---|---:|
| Mean depth | 45.45 |
| Std depth | 20.83 |
| Q1 (25th %) | 29 |
| Q4 (max) | 152 |

Depth distribution is heavy at low depths (Q1 mass dominates by ~10×).
Per-residue Q1 mass dwarfs Q4 mass for all 16 residues. **No residue-
specific Q4-vs-Q1 tilt that correlates with R58 residuals.**

For residues r=5, r=23, r=13 specifically:
- r=5: Q1 mass 3.7M, Q4 mass 241K (ratio 15:1)
- r=23: Q1 mass 2.4M, Q4 mass 247K (ratio 10:1)
- r=13: Q1 mass 2.0M, Q4 mass 243K (ratio 8:1)

The Q1/Q4 ratio doesn't structurally distinguish the QSD-extreme residues.
Esscher tilt parameter λ would need to be different per residue to fix
the residuals — but Esscher tilt is uniform across orbits, so it can't.

## Verdict — outcome (γ): proposed closure REJECTED

R58's residuals at QSD-tilt extremes are **NOT explained by σ-band /
depth-tilt structure**. Esscher tilt by exp(λ · depth) gives at most
+0.014 improvement over baseline Pearson 0.857 — well below the 0.95
target.

### What R58's residuals actually reflect

The baseline residuals follow a pattern consistent with R69's finding
(Chang's P ≠ K kernels): the inverse-tree subtree-size measure and the
forward-trajectory measure (R60 size-stratified Markov, Pearson 0.91)
are **structurally different objects** living on the same support.

R58 averages mass in the inverse tree from m=1; R60 weights by survivor-
conditioned forward orbit visits. These differ by:
1. Backward (R58) vs forward (R60) direction
2. Subtree-size weighting (R58) vs survivor-conditioned visit weighting (R60)
3. Single-resolution residue (R58) vs joint (residue, log-size) (R60)

R69's analysis at the OPERATOR level showed these are not algebraically
related (kernel Frobenius diff 58%, stationary Pearson 0.54). R58's 0.857
Pearson with D_emp is the best the inverse-tree subtree-size measure can
achieve; the gap to D_emp is structural, not σ-tilt-correctable.

### Audit hygiene update

**Claim 6** (Audit / Validation Task 2): "Esscher tilt closes R58's
QSD-tilt residuals" — was classified as "propose-not-demonstrate."

**Update:** **REJECTED.** Esscher tilt at depth (= σ_orbit) does not
improve R58 (+0.014 Pearson at best, far below 0.95 target). The
proposed closure mechanism is wrong.

R58 stays at Pearson 0.857. R60 size-stratified Markov (Pearson 0.91
in-sample, 0.75 train-test) remains the stronger framing of the
trajectory measure. The two frameworks (R58, R60) are complementary
proxies; neither is the "true" trajectory measure (which lives at
infinite resolution and isn't directly computable).

### What this opens

1. **Identify the actual closure path for R58 residuals.** Not σ-tilt;
   maybe joint (residue, depth-bin) marginalization analogous to R60's
   approach but on the inverse tree. Test in follow-up.
2. **R69 framework correspondence:** the R58 vs D_emp gap should be
   explained by the same QSD-vs-invariant-measure distinction R69
   identified for K vs P operators.
3. **Walk back any prior writing** that claims "Esscher tilt closes R58"
   — including in v3.7 framework synthesis chapter if it appears there.

### For v3.7 / external correspondence

> R58's inverse-tree subtree-size measure achieves Pearson 0.857 with
> the trajectory measure D_emp. The largest residuals at r=5, r=23,
> r=13 are STRUCTURAL — they reflect the difference between inverse-tree
> backward weighting and forward-trajectory survivor-conditioned
> weighting (consistent with R69's Chang ≠ K finding at the operator
> level). Tested correction by Esscher tilt at exp(λ · σ_orbit):
> rejected, best λ ≈ 0 with improvement only +0.014.
>
> The R58 framing is a useful proxy for the trajectory measure but
> structurally distinct from it. The R60 size-stratified Markov
> framing (Pearson 0.91 in-sample) is closer; both have explained
> structural residuals at the QSD extremes.

## Files

- `experiments/90_esscher_tilt_r58_closure.py` — full analysis (~3s compute)
- `experiments_output/90_r58_residuals.csv` — per-residue residuals (16 rows)
- `experiments_output/90_esscher_lambda_sweep.csv` — λ sweep (~80 rows)
- `experiments_output/90_sigma_quartile_d_avg.csv` — per-residue Q1..Q4 mass
- `experiments_output/90_esscher_log.txt` — full diagnostic log

## Honest update to STATE.md

Audit Claim 6 should change from:
- ❓ "Esscher tilt closes R58 residuals (propose-not-demonstrate)"

To:
- ❌ "Esscher tilt at σ_orbit does NOT close R58 residuals (R77 demonstrated;
     improvement only +0.014). R58 residuals are structural, reflecting
     inverse-tree vs forward-trajectory measure distinction (R69-style)."
