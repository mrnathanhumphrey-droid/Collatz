# Result 78 — W_visit derivation from R66 3-adic Bohr: outcome (γ); marginal works, per-cell conditional needs additional structure

**Date:** 2026-05-03. R77 follow-up: tests whether the empirical visit-frequency
profile W_visit within each (r mod 32, b) joint state is derivable from R66's
3-adic Bohr stationary π_4 on (Z/81Z)*.

**Verdict (γ):** π_4 explains the **marginal** m mod 81 distribution
to high precision (Pearson 0.987) but NOT the **per-cell conditional**
distribution that R60's empirical kernel actually uses. K_full built with
π_4-derived W_visit gives uniform Perron eigvec — fails to recover D_avg.

## Setup

Walked 1.5M Collatz orbits at N=2^32 (max_T = 250 Syracuse steps), capturing
(r mod 32, b, m mod 81) at each step. Aggregated to compute:
- Empirical W_visit[(r, b)][m mod 81] per joint state
- R66 stationary π_4 on (Z/81Z)* — 54 coprime-to-3 residues
- K_full = K_dynamics weighted by π_4 across the joint (m mod 32, m mod 81)
  → (v, r') transition table

## Step 1 finding: marginal m mod 81 ≈ π_4 (Pearson 0.987)

The trajectory measure's marginal mod 81 is well-described by R66's chain:

```
  Top 10 coprime residues by EMPIRICAL mod-81 mass:
   r81    p_emp     pi_4    ratio
   80    0.0881   0.0897   1.018
   71    0.0711   0.0673   0.947
   26    0.0678   0.0667   0.983
   20    0.0661   0.0668   1.012
   40    0.0460   0.0449   0.975
   ...
```

Pearson(empirical marginal mod 81, π_4) = 0.987 across 54 coprime residues.
MAE = 0.002. The R66 chain stationary IS the trajectory measure's mod-81
marginal.

Mass at m ≡ 0 mod 3 residues: 0.4% (small but nonzero, consistent with
R66's chain ignoring this and accounting for ~99.6% of mass).

## Step 2 finding: per-cell conditional distributions vary widely

Per-(r mod 32, b) cell Pearson with π_4:
- mean = 0.778
- median = 0.895
- min = −0.160 (some cells anti-correlated!)

Sample cells:
```
   r    b    visits   coprime_frac   ρ vs π_4
   3   10   158,593      1.000        0.031   ← FAILS
   5    8   765,336      1.000        0.044   ← FAILS
   1   12   205,613      1.000        0.784   ← partial
  15   20   218,766      1.000        0.998   ← MATCHES
  21    8   115,128      1.000       −0.009   ← FAILS (singular boundary)
  13   10   209,241      1.000        0.407   ← partial
```

**The per-cell conditional distribution is NOT just π_4.** Some cells
(like r=15 b=20) have CRT-independent within-cell distribution that
matches π_4 perfectly. Other cells (r=3 b=10, r=5 b=8, r=21 b=8) have
within-cell concentrations on specific m mod 81 values that depend on
the (r, b) coordinates.

## Step 3 finding: K_full with π_4-derived W_visit gives uniform Perron

```
   r       rho_pred    D_pred     D_avg      diff
   1       0.062541    0.985      1.609     −0.624
   3       0.062411    1.000      1.236     −0.236
   5       0.062437    0.991      1.864     −0.873
   ...    (all rho_pred ≈ 1/16 = 0.0625)
  31       0.062598    1.013      0.767     +0.246

  Total |D_pred − D_avg| = 5.41
  Pearson ρ = −0.71

  R60 empirical reference: total_dev = 3.40, ρ = +0.80
  R77 K_dynamics-only:     total_dev = 5.40, ρ = −0.57
```

K_full ≈ K_dynamics — no improvement. The reason: π_4 is a **global**
weighting. K_full[(r, b) → (r', b')] uses the same π_4 weighting over
m mod 81 for every (r, b) cell. So K_full has translation symmetry in
b within each r, giving uniform Perron after marginalization.

To break this symmetry, W_visit needs (r, b)-conditional weights that
DEPEND on the joint state, not just a global π_4.

## What this means

The trajectory measure's W_visit has structure beyond pure 3-adic Bohr
concentration:

1. **Marginal level (mod 81 alone): π_4 captures it** — R66's chain is
   the right marginal model.
2. **Conditional level (mod 81 given (r mod 32, b)): NOT just π_4** —
   joint (mod 32 × mod 81) coupling matters in a way that breaks CRT
   independence.

The 2-adic and 3-adic structures interact in a (r, b)-dependent way that
isn't captured by either alone. Specifically:
- (r=15 b=20) has CRT-independent structure: m mod 81 is approximately
  uniform within this cell (matching π_4 marginal)
- (r=5 b=8) has strong concentration: m mod 81 is preferential to
  specific values that don't match π_4 ordering

The (r, b)-dependent concentration likely reflects:
- Which higher-bit residues mod 64, 128, etc. are accessible from this
  (r, b) cell under reverse-Syracuse dynamics
- Which 2^k powers fit within the size bin b

A first-principles derivation needs to track the JOINT 2-3-adic
structure of m, not just project to mod 81.

## Per brief outcomes

| Outcome | Status |
|---|---|
| (α) W_visit derived from π_4 | **REJECTED** |
| (β) Partial — some cells derived | **APPLIES** at marginal level + some cells |
| (γ) π_4 alone doesn't explain W_visit | **PRIMARY** |

## Strengthens

- **R65/R66 framework is correct at the marginal level** (Pearson 0.987
  with π_4 marginal). The trajectory measure on Z₂ projects to the
  R66 stationary on (Z/81Z)*.
- **R77's "K_emp = K_dynamics × W_visit" decomposition holds** — the
  W_visit piece is real and not derivable from purely first-principles
  3-adic structure alone.

## Walks back

- **"3-adic Bohr concentration explains W_visit"** (R77 hypothesis)
  is partially true (marginal yes, conditional no). The full W_visit
  structure needs joint 2-3-adic treatment.

## What would close (α)

A first-principles derivation of W_visit needs:
1. Joint (m mod 2^j × m mod 3^k) for j ≥ 6 (beyond mod 32) and k ≥ 4
2. The coupled stationary on this joint state space
3. Project to (r mod 32, b, m mod 81) and check W_visit_predicted

This is a 2592 × 2^j stationary computation — tractable for j up to ~10
(state space ~ 2.6M); larger j needs sparse methods.

Alternative: derive W_visit's b-dependence analytically from log-size
random walk's mixing time across the 3-adic structure.

## Connection to R74 c = 7/45

Both R74 (Fourier-side) and R77/R78 (Markov-side) reduce to:
- Compute leading mode of a deviation-trajectory-measure system
- The trajectory measure's structure couples 2-adic and 3-adic dimensions

**Refined open piece:** the trajectory measure has a JOINT 2-3-adic
Bohr-set structure (not just 3-adic). Pinning that down would close
both R74 (gives c = 7/45) and R77 (gives W_visit) in one analytical
breakthrough.

## Files

- `w_visit_derivation.py` — script
- `w_visit_derivation_log.txt` — full output
- `K_full.npz` — sparse K_full (with π_4 weighting)
- `K_full_perron_results.csv` — Perron marginal vs D_avg
- `w_visit_predicted_vs_empirical.csv` — per-cell Pearson with π_4
- `w_visit_derivation.md` — this writeup

## Concrete next moves

1. **Joint 2-3-adic stationary**: build Markov chain on (Z/64Z × Z/81Z)*
   joint state and compute stationary; compare W_visit projection to
   empirical at deeper modulus.
2. **(r, b)-conditional W_visit empirical fit**: extract empirical
   W_visit[(r, b)][m mod 81] per cell and use IT in K_full → check
   Perron recovers D_avg. This isolates whether the per-cell
   conditional is the load-bearing piece (vs the marginal).
3. **Cell-class identification**: cluster (r, b) cells by their
   per-cell Pearson with π_4 — which cells have CRT-independent
   structure (Pearson > 0.95) vs concentrated structure (Pearson
   < 0.5)? Map to size-bin / residue patterns.
4. **Update STATE.md** open piece 9 with this refinement.
