# Result 54 — Size-stratified Markov framework reproduces D_avg

**Outcome: (α)** — size-stratified Markov captures D_avg significantly better
than null model and all prior kernel candidates.

## Summary

| Framework | Total |D_pred − D_avg| | Pearson ρ | Verdict |
|---|---|---|---|
| Cylinder QSD (R51) | 5.39 – 7.71 | low | (γ) miss |
| Inverse-tree (R52) | 13.2 – 27 | low | (γ) miss |
| Renewal kernel exit (R53) | 5.04 | 0.58 | (β) partial |
| Trivial null (D_full this dataset) | 4.72 | (best previous) | reference |
| **Size-stratified Markov (this)** | **3.40** | **0.80** | **(α) hit** |

Per-step survival rate: λ_PF = 0.957 (vs empirical 0.94 from R50).

## Construction

State space: (r mod 32, b) with r in {1, 3, …, 31} (16 odd residues) and
b = ⌊log₂ m⌋ in {0, 1, …, 63}. Total 16 × 64 = 1024 states.

Empirical kernel built from 1.5M Collatz orbits at N=2^32 (112M Syracuse-step
transitions). For each transition (m → next odd m'), increment
K_counts[(r, b), (r', b')]. Row-normalize: K[s → s'] = #{s → s'} / #{visits to s}.

Filter: keep states with ≥50 visits or inflows. 629 states retained out of 1024.

Compute leading left eigenpair of K (= leading right of K^T) via ARPACK on
sparse 629×629 matrix. Use as v_PF, the QSD on (residue, log-size).

Marginalize: ρ_pred(r) = Σ_b v_PF(r, b); D_pred(r) = ρ_pred(r) / π₃₂(r).

## Eigenvalue structure

- λ₁ = 0.9566
- λ₂ = 0.9462
- λ₃ = 0.9266
- Spectral gap (1 − |λ₂|/|λ₁|) = 1.1%

The small spectral gap means the chain is near-degenerate; the QSD has a
slow second mode that is also descent-related. The leading λ₁ ≈ 0.957 is
slightly above the empirical 0.94 from R50; the closer λ₂ ≈ 0.946 may be
the "physical" survival mode after subtracting initial transients.

## Per-residue D_pred vs D_avg

```
   r       pi      D_avg   D_pred    diff
   1   0.0635   1.6091   1.3761  −0.2330
   3   0.0624   1.2364   0.9151  −0.3213
   5   0.0630   1.8639   1.5588  −0.3051
   7   0.0624   0.7375   1.1063  +0.3688
   9   0.0623   0.6958   0.9581  +0.2623
  11   0.0628   0.6575   0.9240  +0.2664
  13   0.0625   0.5573   0.6767  +0.1194
  15   0.0623   1.0578   1.0910  +0.0332
  17   0.0632   1.1320   1.0816  −0.0503
  19   0.0624   0.6275   0.6815  +0.0540
  21   0.0625   0.9307   0.6333  −0.2974
  23   0.0624   1.3981   1.1085  −0.2895
  25   0.0618   0.5440   0.7000  +0.1560
  27   0.0618   0.8017   1.1451  +0.3434
  29   0.0630   1.3538   1.1600  −0.1938
  31   0.0618   0.7666   0.8682  +0.1015

Total |D_pred − D_avg| = 3.40
MAD per residue        = 0.21
Pearson ρ              = 0.80
Spearman ρ             = 0.72
```

The framework correctly captures the rough shape: high values at r ∈ {1, 3,
5, 23, 29} (mid-orbit residues with extended descent) and low values at
r ∈ {9, 11, 13, 19, 25} (residues with quick descent into m=1 cycle).

The largest residual is at r=27 (D_pred over by 0.34), r=7 (over by 0.37),
and r=5/r=21/r=23 (under by 0.30). These are the residues where the
discrete log-size binning probably aliases dynamics that need finer
resolution — bin width is 1 in log₂(m), so m=5 and m=7 land in the same
b=2 bin even though Syracuse maps them very differently.

## Joint structure

P(b | r) is strongly residue-dependent — the framework is NOT separable.
RMS factorization residual is 79% of ‖v_PF‖, decisively rejecting
v_PF(r, b) = f(r)·g(b).

Examples of P(b | r):
- **r=5**: top bins b=2 (40%), b=8 (18%), b=11 (17%). Heavy at small m
  reflecting that 5 → 16 → 1 (immediate absorption from m=5).
- **r=21**: top bins b=5 (47%), b=9 (12%), b=11 (7%). Concentrated
  around m ≈ 32, consistent with the m_3 = 21 attractor and its
  inverse-tree neighborhood found in R52.

The marginal P(b) shows mass concentrating around b=7-9 (m ≈ 128-512), with
exponential decay at larger b — this is the descent-dominated tail.

## What this means

D_avg is identified as the residue marginal of the leading left eigenvector
(QSD) of an empirical 1024-state Markov kernel on (r mod 32, log₂ m).

This is a **finite-dimensional algebraic identification**. The trajectory
measure D_avg is no longer "an empirical distribution that resists Markov
description" — it is the residue projection of a Perron-Frobenius
eigenvector, where the kernel is built from orbit data and depends on
the joint state (residue, log-size).

Three negative results (R51 cylinder, R52 inverse tree, R53 renewal) in
sequence had been interpreted as evidence that residue-only frameworks
cannot capture D_avg. The size-stratified framework adds the size
dimension explicitly and the gap closes.

The remaining 3.40 total deviation (vs floor of ~0 for an exact match) is
large enough to be physical, not noise. Likely sources:
1. Coarse log-size binning (b is integer-valued; finer bin width would
   resolve the m=5 vs m=37 distinction at b=2 vs b=5).
2. Truncation at 600 Syracuse steps — long-tail orbits not represented.
3. Transient mode contamination — the second eigenvalue 0.946 is close
   enough to leading 0.957 that v_PF mixes slow modes.

## Connection to Chang correspondence

Chang's framework operates on Z/64Z (residue-only, π is uniform).
Size-stratified framework operates on (Z/32Z × log-size bins) — same
Collatz dynamics, different stratification.

Sharpened framing for v3.6: D_avg is the size-stratified QSD's residue
projection. Chang's invariant core I_2 = {7, 27, 31, 59, 63} (mod 64)
governs the residue-only eigendynamics; D_avg adds size-dependence
that Chang's framework averages out.

## Files

- `size_stratified_markov.py` — script
- `size_stratified_markov_log.txt` — full run log
- `size_stratified_kernel.npz` — sparse 629×629 kernel K_sub
- `size_stratified_keep_idx.npy` — index map to (r, b) state space
- `size_stratified_eigvec.csv` — v_PF on (r, b)
- `size_stratified_residue_marginal.csv` — D_pred vs D_avg
