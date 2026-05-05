# R58 closure attempt #2 — Esscher tilt at log(R58/R60) — VERDICT: γ (closure FAILS)

**Date:** 2026-05-05. Companion to [esscher_tilt_r58_closure.md](esscher_tilt_r58_closure.md) (the prior σ_orbit attempt, also rejected). Tests R69's mechanism: the per-node log-ratio of R58 vs R60 weightings as the closure observable.

## Verdict

> **γ — closure FAILS to reach 0.95 target.** The geometric interpolation w_λ(m) = w_R58(m)^(1−λ) · w_R60(m)^λ peaks at **λ=0.95 with Pearson +0.8666** (improvement only +0.010 over R58 baseline 0.8568), with **MAE 0.114** (vs baseline 0.119). Train-test at N=2^21 → N=2^22 gives identical λ_best=0.95 and identical Pearson — no overfitting in either direction.
>
> **Per-residue: tilt fixes r=5/13/23 BUT breaks r=1 and r=21.** The QSD-extremes improve substantially (r=5: +0.46→+0.24, r=13: −0.21→−0.04, r=23: +0.31→+0.27) but two previously-clean residues (r=1: −0.02→−0.23, r=21: −0.002→+0.19) deteriorate by comparable magnitude. The improvements come from REDISTRIBUTING residual mass, not eliminating it.
>
> **R69's mechanism is partially supported but insufficient.** The per-node weight ratio IS structurally meaningful (the geometric interpolation has a clear peak), but R58 and R60 share enough common structural gap with D_emp that mixing them within their convex hull cannot reach 0.95. The R58/R60 gap is NOT the full obstruction.

## Setup

For each odd node m in the inverse tree from m=1 (built at N=2^22, 1.247M nodes):

- **w_R58(m) = subtree_size(m)** — count of descendants, used by R58's inverse-tree subtree-size measure (Pearson 0.857)
- **w_R60(m) = π_R60(r, b) / N_cell(r, b)** — R60's size-stratified Markov stationary mass on cell (m mod 32, ⌊log₂ m⌋), normalized by cell occupancy so that summing over a cell recovers π_R60(r, b). This makes λ=1 reproduce R60's residue-marginal prediction (modulo cell-occupancy weighting).

The closure observable is the per-node log-ratio r(m) = log(w_R58(m) / w_R60(m)). Esscher tilt by r(m) is equivalent to geometric interpolation:

  w_λ(m) = w_R58(m) · exp(λ · r(m)) = w_R58(m)^{1−λ} · w_R60(m)^λ

Computed in log-space with max-subtraction for numerical stability.

## Stage 1+2: Baseline reproduction

Reproduces both R58 (Pearson 0.8568) and R60 (Pearson 0.8574) at N=2^22.

**Note:** R60's "Pearson 0.91" reported in `size_stratified_markov.md` is the joint (r, b) correlation, NOT the residue-marginal one. At the residue-marginal level R60 = R58 to within 0.001 — they're essentially equivalent in this restricted observable.

## Stage 3: Per-residue log-ratio distribution

| r | n | mean log-ratio | std | skew | D_emp − D_R58 |
|---:|---:|---:|---:|---:|---:|
| 1  | 100,620 | +20.24 | 1.40 | −0.87 | −0.02 |
| 3  | 83,132  | +20.02 | 1.46 | −0.64 | +0.02 |
| **5**  | **121,281** | **+20.49** | 1.38 | −1.38 | **+0.46 ← QSD-extreme** |
| 7  | 40,035  | +19.37 | 1.68 | −0.23 | −0.10 |
| 9  | 46,203  | +19.36 | 1.63 | −0.17 | −0.12 |
| 11 | 73,950  | +19.89 | 1.47 | −1.06 | −0.13 |
| **13** | **121,264** | **+20.50** | 1.39 | −0.97 | **−0.21 ← QSD-extreme** |
| 15 | 30,805  | +19.13 | 1.82 | −0.02 | +0.10 |
| 17 | 121,295 | +20.50 | 1.38 | −1.36 | −0.14 |
| 19 | 73,957  | +19.88 | 1.48 | −0.51 | −0.05 |
| 21 | 128,117 | +20.57 | 1.39 | −1.01 | −0.00 |
| **23** | **52,358** | **+19.72** | 1.63 | −0.98 | **+0.31 ← QSD-extreme** |
| 25 | 92,401  | +20.14 | 1.41 | −0.75 | −0.14 |
| 27 | 46,201  | +19.37 | 1.63 | −0.18 | −0.04 |
| 29 | 100,664 | +20.23 | 1.39 | −1.33 | −0.01 |
| 31 | 15,423  | +18.52 | 2.03 | +0.36 | +0.05 |

**Critical observation: the log-ratio means do NOT track the residuals.** r=5 and r=13 have NEAR-IDENTICAL mean log-ratios (+20.49 and +20.50), but R58's residuals at these residues are **opposite-sign** (+0.46 vs −0.21). r=21 has the HIGHEST mean log-ratio (+20.57) but a residual of essentially zero. R69's hypothesis "log-ratio explains the residuals via uniform tilt" is **falsified at the per-residue level** — the observable doesn't carry the right information to distinguish enhancement from depletion.

## Stage 4: Lambda sweep

Pearson on geometric-interpolation w_λ across λ ∈ [0, 1] step 0.05:

| λ | Pearson | MAE |
|---:|---:|---:|
| 0.00 (R58) | +0.8568 | 0.1189 |
| 0.10 | +0.8483 | 0.1202 |
| 0.20 | +0.8352 | 0.1244 |
| 0.30 | +0.8243 | 0.1288 |
| 0.40 | +0.8182 | 0.1318 |
| **0.45** | **+0.8170 (min)** | 0.1327 |
| 0.50 | +0.8172 | 0.1331 |
| 0.60 | +0.8214 | 0.1323 |
| 0.70 | +0.8308 | 0.1294 |
| 0.80 | +0.8457 | 0.1241 |
| 0.90 | +0.8628 | 0.1172 |
| **0.95 (max)** | **+0.8666** | **0.1142** |
| 1.00 (R60) | +0.8574 | 0.1168 |

The Pearson surface is **U-shaped over [0, 1]**: it drops from 0.857 at λ=0 to a trough of 0.817 at λ=0.45, then climbs back through 0.86 at λ=0.95 and falls slightly to 0.857 at λ=1. The minimum at λ≈0.45 means **R58 and R60 disagree on which residues to enhance/deplete in opposite directions** — mixing them produces worse predictions than either alone over much of the interpolation range.

This is a strong signal that **R58 and R60 are not nested or convex-related proxies for the same target**. They each capture some structure but their disagreement is structural, not refinable by interpolation.

## Stage 5: Per-residue residuals at λ_best = 0.95

| r | D_emp | D_R58 | D_R60 | D_λ=0.95 | res_R58 | res_best | |Δ|res| |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.217 | 1.240 | 1.400 | 1.444 | −0.022 | **−0.226** | **+0.20 worse** |
| 3 | 0.990 | 0.972 | 0.914 | 0.919 | +0.018 | +0.071 | +0.05 worse |
| **5** | 1.836 | 1.376 | 1.573 | **1.598** | **+0.460** | **+0.239** | **−0.22 better** |
| 7 | 0.912 | 1.014 | 1.104 | 1.107 | −0.101 | −0.194 | +0.09 worse |
| 9 | 0.794 | 0.915 | 0.955 | 0.922 | −0.121 | −0.127 | ~ |
| 11 | 0.849 | 0.979 | 0.928 | 0.929 | −0.129 | −0.079 | −0.05 better |
| **13** | 0.680 | 0.889 | 0.676 | **0.718** | **−0.209** | **−0.038** | **−0.17 better** |
| 15 | 1.093 | 0.989 | 1.087 | 1.047 | +0.104 | +0.047 | −0.06 better |
| 17 | 1.010 | 1.151 | 1.095 | 1.099 | −0.140 | −0.089 | −0.05 better |
| 19 | 0.771 | 0.816 | 0.680 | 0.680 | −0.045 | +0.091 | +0.05 worse |
| 21 | 0.862 | 0.864 | 0.632 | 0.669 | −0.002 | **+0.193** | **+0.19 worse** |
| **23** | 1.351 | 1.043 | 1.107 | **1.086** | **+0.308** | **+0.265** | −0.04 better |
| 25 | 0.702 | 0.842 | 0.691 | 0.706 | −0.140 | −0.003 | −0.14 better |
| 27 | 0.967 | 1.005 | 1.133 | 1.070 | −0.038 | −0.103 | +0.07 worse |
| 29 | 1.112 | 1.124 | 1.170 | 1.171 | −0.011 | −0.059 | +0.05 worse |
| 31 | 0.835 | 0.783 | 0.855 | 0.836 | +0.052 | −0.001 | −0.05 better |

**The tilt redistributes residual mass — it doesn't reduce total error magnitude.** Sum of absolute residuals:
- R58 (λ=0): **1.90**
- λ=0.95: **1.83** (4% improvement)

Improvements concentrated at r=5/13/25 (cumulative −0.53 absolute reduction); deteriorations concentrated at r=1/21/7 (cumulative +0.49 absolute increase). Net: marginal.

QSD-extreme verdict: r=5 and r=13 ARE materially improved (residuals halved or better). r=23 essentially unchanged. **r=1 and r=21 are NEW QSD-extreme failures** induced by the tilt — they had residuals near zero in R58 and now sit at ±0.2.

## Stage 6: Train-test (no overfitting)

| metric | value |
|---|---:|
| λ_best at N=2^21 (train, 622,844 nodes) | 0.95 |
| Pearson_train | +0.8665 |
| Apply λ_train at N=2^22: Pearson_test | +0.8666 |
| Train-test gap | 0.0001 |

Identical λ_best across two N values. The structural relationship is N-stable; the tilt isn't capturing N-specific noise. But it's not capturing the underlying mechanism either.

## Why the closure fails

Three structural reasons emerge from the data:

1. **Per-node log-ratio doesn't carry sign information about which way to push each residue class.** r=5 (over-predicted depletion, residual +0.46) and r=13 (over-predicted enhancement, residual −0.21) have nearly identical mean log-ratios (+20.49 vs +20.50). A uniform Esscher tilt at log(R58/R60) cannot distinguish them. Same as the σ_orbit attempt's failure mode, just at a different observable.

2. **R58 and R60 are equivalent proxies at residue-marginal resolution** (Pearson 0.8568 vs 0.8574). The 0.91 advantage R60 claimed elsewhere is at the joint (r, b) level. At the marginal level, R60 doesn't beat R58 enough to give the convex combination room to outperform either endpoint.

3. **The interpolation U-shape (min Pearson 0.817 at λ=0.45) shows R58 and R60 disagree in OPPOSITE directions on residual signs.** Their gap with D_emp isn't a magnitude difference (which mixing could fix); it's a sign-pattern difference (which mixing makes worse before getting back to baseline at the endpoints).

The interpretation is **R69's mechanism is real but operationally insufficient**: the weight ratio captures part of the structure (improving r=5/13/23 substantially) but not the full pattern (breaking r=1/21). The full closure mechanism requires per-residue-specific tilts, not a single λ.

## Decision-rule outcome

Per the brief's decision tree:

- ❌ **(A) λ_best ∈ [0.4, 1.0] with Pearson > 0.95 AND r=5/13/23 simultaneously improved**: Pearson never reaches 0.95. r=5/13/23 improved but other residues broke.
- ❌ **(B) λ_best near 1 AND Pearson 0.91 (matches R60 baseline)**: Doesn't match R60's joint-level 0.91; matches the residue-marginal R60 of 0.857.
- ❓ **(C) λ_best fails to simultaneously fix r=5/13/23**: r=5/13/23 individually improved, but r=1 and r=21 broke. Mixed outcome.
- ✅ **(D) Pearson never exceeds 0.92**: Confirmed. Best Pearson 0.8666 < 0.92.

**Best classification: (D) with partial-(C) qualifier.** The convex combination of R58 and R60 doesn't reach the closure threshold; both proxies share structural gaps with D_emp that aren't algebraically related; expected from R69's stationary Pearson 0.54 at the operator level (not algebraically related kernels can't be linearly combined to recover the true measure).

## What this opens (and what it closes)

**Closes:**
- Both Esscher-tilt closure attempts (σ_orbit and log(R58/R60)) for R58's QSD-extreme residuals: REJECTED.
- Hypothesis "R58/R60 weight ratio is the full obstruction": **falsified**. The gap has additional structure beyond what either proxy captures.

**Opens — what the closure path actually requires:**

The two failed Esscher attempts (σ_orbit and log(R58/R60)) share a single failure mode: a uniform tilt observable cannot carry per-residue sign information that residue-marginal residuals require. r=5 (residual +0.46) and r=13 (residual −0.21) need opposite-direction corrections; any uniform λ pushes them in the same direction. The structural takeaway is that **closure requires non-uniform tilt**, with three sub-options:

- **Non-uniform tilt** — residue-conditional λ_r per residue. 16 free parameters on 16 datapoints is empty fitting; not a real test, but worth recording as the pattern.
- **A different observable** that varies sign-aligned with residuals across residues. Neither σ_orbit nor log(R58/R60) does. Identifying such an observable analytically is open.
- **A different closure mechanism entirely** — outside the Esscher-tilt family. R58 stays at Pearson 0.857 as the best the inverse-tree subtree-size measure can achieve; the closure mechanism may not pass through tilting at all.

**Erratum on a prior framing.** An earlier draft of this writeup cited the joint 2-3-adic Bohr empirical positive (`result_bohr_probe.md`, z=16.5 at k=20) as the "natural successor probe." That citation is **withdrawn**: the Bohr aggregate signal was deflated 2026-05-04/05 by the bracket-stratification probe ([result_bohr_probe_strat.md](result_bohr_probe_strat.md)). Per-bracket chi²/df at v ∈ (10⁶, 10⁹] is 0.95 (z = −0.99) and at v > 10⁹ is 0.94 (z = −1.18) — statistically CRT-independent within ±2σ. The original aggregate signal was driven by the v ≤ 100 bracket-A descent funnel (low-v trajectories transiting toward 1), not joint structure at scales relevant to D_emp. The "3-adic structure that residue-mod-32 erases" hypothesis is not currently supported by an empirical observable.

## Files

- [esscher_tilt_r58_closure_v2_log_ratio.py](esscher_tilt_r58_closure_v2_log_ratio.py) — script
- [esscher_tilt_v2_log_ratio_lambda_sweep.csv](esscher_tilt_v2_log_ratio_lambda_sweep.csv) — full λ sweep, 21 rows
- [esscher_tilt_v2_log_ratio_residuals.csv](esscher_tilt_v2_log_ratio_residuals.csv) — per-residue at λ=0 vs λ=0.95
- [esscher_tilt_v2_log_ratio_perclass_logratio.csv](esscher_tilt_v2_log_ratio_perclass_logratio.csv) — log-ratio distribution stats per residue
- [esscher_tilt_r58_closure_v2_verdict.md](esscher_tilt_r58_closure_v2_verdict.md) — this writeup

## STATE.md impact

Update Audit Claim 6 to:
> **REJECTED at both observables tested.** Esscher tilt at σ_orbit (depth) gives +0.014 Pearson improvement; Esscher tilt at log(R58/R60) per R69's mechanism gives +0.010 improvement. Neither reaches 0.95. R58 stays at Pearson 0.857 as the best inverse-tree subtree-size measure. The structural takeaway across both attempts: **a uniform Esscher tilt observable cannot carry per-residue sign information when residuals at different residues demand opposite-direction corrections.** Closure requires non-uniform tilt, a sign-aligned observable, or a different mechanism entirely. R69's weight-ratio mechanism is partial and operationally insufficient.
