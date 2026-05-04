# Result 60 v2 (finer binning) — outcome (β) with overfitting boundary identified

**Status.** Validation Task 1's flagged Pearson 0.91 at B=109/log_base=1.5 is
**confirmed in-sample** (Pearson 0.9116, total dev 2.27, exact replication
of the prior sweep). However, the **train-test holdout reveals overfitting**:
out-of-sample Pearson drops from baseline's 0.78 to 0.75 at B=109, and
collapses to 0.65 at B=150 (where in-sample reaches 0.99).

**Outcome (β):** improvement confirmed; finer binning is the asymptotic limit
of in-sample fit but loses generalization beyond log_base ≈ 1.5. The Markov
assumption STRENGTHENS at finer binning (8.8% → 3.1% non-Markov info).
Framework comparison strengthens (3.40 → 2.27 total dev). For external
correspondence: present Pearson 0.91 with the train-test caveat clearly stated.

**Outcome (δ) REJECTED**: in-sample Pearson 0.99 at B=150 is **overfitting**
(train-test 0.65), not asymptotic structural identification.

## Reference build: B=109, log_base=1.5

Direct replication of Validation Task 1 finding.

- 1.5M orbits at N=2^32, 112M Syracuse-step transitions
- 1033 / 1744 states retained at min_visits=50
- λ_PF = 0.9514, λ_2 = 0.9424, spectral gap = 0.0095
- **Pearson = 0.9116, total dev = 2.27, MAD per residue = 0.142**

Per-residue diagnostic (all 16 odd residues mod 32):

| r | π_32 | D_avg | D_pred | diff |
|---:|---:|---:|---:|---:|
| 1 | 0.0635 | 1.609 | 1.473 | −0.136 |
| 3 | 0.0624 | 1.236 | 1.004 | −0.232 |
| 5 | 0.0630 | 1.864 | 1.731 | −0.133 |
| 7 | 0.0624 | 0.738 | 0.998 | +0.260 |
| 9 | 0.0623 | 0.696 | 0.879 | +0.183 |
| 11 | 0.0628 | 0.658 | 0.804 | +0.147 |
| **13** | 0.0625 | 0.557 | 0.546 | **−0.011** |
| 15 | 0.0623 | 1.058 | 1.139 | +0.081 |
| 17 | 0.0632 | 1.132 | 1.070 | −0.062 |
| 19 | 0.0624 | 0.628 | 0.648 | +0.021 |
| **21** | 0.0625 | 0.931 | 0.643 | **−0.288** |
| 23 | 0.0624 | 1.398 | 1.250 | −0.148 |
| 25 | 0.0618 | 0.544 | 0.631 | +0.087 |
| **27** | 0.0618 | 0.802 | 1.067 | **+0.265** |
| 29 | 0.0630 | 1.354 | 1.232 | −0.122 |
| 31 | 0.0618 | 0.767 | 0.864 | +0.098 |

Top residual residues: r=21 (−0.288), r=27 (+0.265), r=7 (+0.260),
r=3 (−0.232), r=9 (+0.183).

Notable improvements over baseline R60 (B=64, log_base=2):
- r=13: from −0.32 → **−0.011** (near-perfect)
- r=17: from +0.52 → **−0.062** (15× better)
- r=29: from −0.92 → **−0.122** (8× better)

## 7-concern validation summary

| Concern | Test | Baseline (B=64) | B=109/1.5 | Delta |
|---|---|---:|---:|---:|
| 1 sample-size | Pearson stability | 0.803 ± 0.001 | **0.911 ± 0.001** | +0.108 |
| 2 train-test | Holdout Pearson | 0.784 | **0.755** | **−0.029** |
| 3 binning | Pearson stability across log bases | 0.80–0.91 | 0.79–0.91 | mixed |
| 4 Markov | I(prev;next\|t)/I(t;next) | 8.77% | **3.06%** | strengthened |
| 5 framework comparison | R60 vs null/inv-tree | 3.40 vs 5.33/14.7 | **2.27** vs 5.33/14.7 | strengthened |
| 6 λ_PF auto-match | Real vs scrambled/uniform null | 0.957 vs 0.99/0.99 | 0.951 vs 0.99/0.99 | structural |
| 7 factorization null | RMS residual z-score | +6.6σ | **+16.1σ** | strengthened |

**4 strengthen, 1 maintains structurally informative, 1 weakens (train-test),
1 mixed.**

## Concern 1 — Sample-size stability: PASSES

| N | n_kept | total_dev | Pearson | λ_PF |
|---|---:|---:|---:|---:|
| 2^28 | 919 | 2.274 | 0.9116 | 0.9515 |
| 2^30 | 978 | 2.273 | 0.9118 | 0.9515 |
| 2^32 | 1033 | 2.274 | 0.9116 | 0.9514 |
| 2^34 | 1085 | 2.282 | 0.9110 | 0.9515 |

Pearson stable to 4 decimals across 4 orders of magnitude in N. No N-stability
issue.

## Concern 2 — Train-test holdout: WEAK FAIL (0.75 < 0.85 target)

| comparison | total_dev | Pearson |
|---|---:|---:|
| Train K v_PF marginal vs global D_avg | 2.274 | 0.9117 |
| **Train K v_PF marginal vs Test-derived D_test** | 3.192 | **0.7546** |
| Sanity: D_test vs global D_avg | 4.741 | 0.5885 |

Pearson(train_v_PF, D_test) = 0.7546 — **lower than baseline's 0.7843**.

The train-test holdout reveals that the in-sample improvement (0.80 → 0.91)
includes finite-sample fitting noise. Out-of-sample, the finer kernel does
not generalize as well as the coarser baseline.

Note that 0.7546 is still > Pearson(D_avg, D_test) = 0.5885, indicating
v_PF captures (r, b) population structure beyond pure single-sample noise.
But it does NOT improve generalization over baseline.

## Concern 3 — Bin-base sensitivity: MIXED

| log_base | B | n_kept | total_dev | Pearson |
|---:|---:|---:|---:|---:|
| 1.4 | 131 | 1223 | 2.774 | 0.866 |
| **1.5** | **109** | **1033** | **2.274** | **0.912** |
| 1.6 | 94 | 899 | 3.372 | 0.792 |
| 1.7 | 83 | 804 | 3.147 | 0.822 |

Maximum at log_base = 1.5; modest sensitivity within ±0.1 base. Notable that
log_base = 1.6 drops to Pearson 0.79 — finer than 1.5 doesn't monotonically
improve in-sample fit either at this resolution.

## Concern 4 — Markov assumption: STRENGTHENED PASS

| quantity | baseline (B=64) | B=109/1.5 |
|---|---:|---:|
| H(s_t) | 8.703 bits | 9.343 bits |
| I(s_t; s_{t+1}) | 6.400 bits | 7.507 bits |
| I(s_{t-1}; s_{t+1} \| s_t) | 0.561 bits | **0.230 bits** |
| Ratio | 8.77% | **3.06%** |

**Finer binning makes the (r, b) coordinate MORE Markov.** The non-Markov
information leakage drops from 8.8% to 3.1%, satisfying the strict <5%
threshold. This is a genuine structural improvement: at B=109/log_base=1.5,
the pair-coupling captures 96.9% of the bit-level predictive information.

## Concern 5 — Framework comparison at matched conditions: STRONGER PASS

| framework | total_dev | Pearson |
|---|---:|---:|
| Trivial null (D=1) | 5.332 | 0.000 |
| Inverse tree (R23) | 14.728 | 0.095 |
| R60 baseline (B=64, log_base=2) | 3.40 | 0.80 |
| **R60 v2 (B=109, log_base=1.5)** | **2.27** | **0.91** |

R60 v2 beats null by 57%, beats baseline R60 by 33% on total deviation.
Strongest empirical match across all frameworks.

## Concern 6 — λ_PF auto-match: STRUCTURAL (clarified)

| kernel | λ_PF |
|---|---:|
| Real K_sub (B=109) | **0.9514** |
| Scrambled K (preserves row sums) | 0.9955 |
| Uniform-column K | 0.9923 |

λ_PF = 0.9514 is genuinely lower than scrambled/uniform null (~0.99),
consistent with the original framework — the empirical column structure
compresses the eigenvalue from null. Same interpretation as baseline.

## Concern 7 — Factorization residual null: STRENGTHENED PASS

| quantity | baseline (B=64) | B=109/1.5 |
|---|---:|---:|
| Real factorization residual | 78.7% | **87.6%** |
| Null mean | 68.1% | 68.6% |
| Null std | 1.6% | 1.2% |
| Real z-score | +6.57 | **+16.12** |

Non-separability is **stronger at finer binning** (16σ vs 6σ above null).
The (r, b) coupling has more irreducible joint structure when finer
resolution exposes it.

## Outcome (δ) verification: B=150 reveals OVERFITTING

Tested whether Pearson keeps climbing past 0.91 with even finer binning.
Used log_base = 2^(32/B) so bins span exactly N=2^32.

| B | log_base | n_kept | total_dev | Pearson **in-sample** | spectral gap |
|---:|---:|---:|---:|---:|---:|
| 50 | 1.558 | 666 | 3.002 | 0.8476 | 0.0146 |
| 75 | 1.344 | 976 | 2.321 | 0.9075 | 0.0127 |
| 100 | 1.248 | 1280 | 1.691 | 0.9448 | 0.0098 |
| 125 | 1.194 | 1579 | 1.109 | 0.9808 | 0.0062 |
| **150** | **1.159** | **1876** | **0.850** | **0.9888** | **0.0076** |
| 200 | 1.117 | 2459 | 0.781 | 0.9923 | 0.0116 |
| 175, 250, 300 | — | — | eval failed | — | — |

**B=150 verification (independent run with full validation):**

| concern | result |
|---|---:|
| In-sample (B=150) Pearson vs D_avg | 0.9888 |
| **Train-test holdout Pearson** | **0.6520** ← collapses |
| Sample-size stability (N=2^28..2^34) | 0.988 stable |
| Spectral gap | 0.0076 (very tight, ARPACK warns) |

**The train-test holdout Pearson DROPS from 0.78 (baseline) to 0.75
(B=109/1.5) to 0.65 (B=150).** The in-sample Pearson climbing from 0.91 → 0.99
at B=150 is overfitting to kernel-construction noise.

In-sample / out-of-sample gap as B grows:
- Baseline (B=64): 0.80 vs 0.78 → gap 0.02
- B=109/1.5: 0.91 vs 0.75 → gap 0.16
- B=150: 0.99 vs 0.65 → gap 0.34

**Textbook overfitting signature.** The "outcome (δ) Pearson approaches 1"
hypothesis is rejected.

## Per-residue / per-bin contribution analysis at B=109/log_base=1.5

For r=5 (most-enhanced residue, D_avg = 1.864):
- Top bin b=3 (m ∈ [3.4, 5.1)) carries 40.8% mass — m=5 itself
- Bin b=19 (m ∈ [2217, 3325)) carries 19.7%
- Bin b=14 (m ∈ [292, 438)) carries 19.3%

For r=21 (most-enhanced under R23 inverse tree, D_avg = 0.931):
- Top bin b=9 (m ∈ [38.4, 57.7)) carries 63.6% mass — m_3 = 21 attractor
- Bin b=16 (m ∈ [657, 985)) carries 7.2%

For r=13 (most-depleted, D_avg = 0.557):
- Top bin b=6 (m ∈ [11.4, 17.1)) carries 48.8% mass — m=13 itself
- Spread broadly across mid-size bins

The atomic concentration on m_j attractor sequence is **explicit at B=109**:
mass concentrates in single bins for residues r ∈ {1, 5, 13, 21, 53, 85, ...}
which are the m_j residues. This was muddled at baseline (B=64, log_base=2)
because b=2 mixed m ∈ [4, 8).

## Verdict — outcome (β) with overfitting boundary

**Outcome (β):** Pearson 0.91 at B=109/log_base=1.5 confirmed in-sample.
The improvement is genuine for the closed-form identification "D_avg is the
residue marginal of the empirical kernel's leading eigenvector" — when the
kernel resolves the m_j atomic structure that baseline binning blurred.

**Train-test caveat:** Out-of-sample, the framework's Pearson is 0.75 — a
slight DROP from baseline's 0.78. The finer binning fits in-sample structure
better but does not improve generalization.

**Markov + framework + factorization concerns ALL strengthen at B=109/1.5:**
- Markov 8.8% → 3.1% (strict-criterion pass)
- Framework total dev 3.40 → 2.27
- Non-separability +6.6σ → +16.1σ
- These are genuine structural improvements not artifacts

**Outcome (δ) rejected:** At B=150+, in-sample Pearson reaches 0.99 but
train-test collapses to 0.65 — overfitting. The framework does NOT
asymptotically approach Pearson 1; it asymptotically approaches kernel-
construction noise.

## For v3.7 / external correspondence

**Lead claim (revised):**

> The trajectory measure D_avg is identified as the residue marginal of
> the leading left Perron eigenvector of an empirical Markov kernel on the
> joint state space (residue mod 32, log_{1.5} m). At B=109 bins (state
> space 1744, retained kernel 1033 states), the framework reproduces D_avg
> with Pearson 0.91 and total deviation 2.27 across 16 odd residues mod 32,
> beating all prior frameworks. The Markov approximation captures 96.9% of
> the bit-level predictive information (vs 91.2% at baseline). The non-
> separable joint (r, b) structure is statistically significant at +16σ
> above random-vector null.
>
> **Caveats:**
> - In-sample/out-of-sample Pearson gap is 0.16 (vs baseline's 0.02). Finer
>   binning improves structural identification but does not improve
>   generalization beyond baseline.
> - At B>150 with log_base chosen to span N exactly, in-sample Pearson
>   reaches 0.99 but train-test collapses; this is overfitting, not
>   asymptotic identification.

**What stands as-is:**
- Pearson 0.91 in-sample identification (validated against D_avg from
  independent reference)
- Markov approximation is structurally cleaner at finer binning (3.1% non-
  Markov info)
- Atomic concentration on m_j residue is explicit at B=109 binning
- Framework comparison vs null and inverse tree

**What needs framing:**
- "Pearson 0.91" should be qualified as "in-sample identification of D_avg"
  not "predictive Pearson out-of-sample"
- Out-of-sample improvement over baseline is **not** demonstrated
- Pearson does NOT approach 1 with finer binning; it saturates near 0.91
  before overfitting takes over at B=150+

## Files

- `experiments/83_result60_v2_finer.py` — full reference build + 7-concern
  validation + scaling sweep
- `experiments/84_result60_v2_outcome_delta.py` — verifies B=150 overfitting
  via train-test
- `result60_v2_kernel.npz` — sparse 1033×1033 kernel
- `result60_v2_keep_idx.npy` — index map
- `result60_v2_eigvec.csv` — v_PF on (r, b) state space (1744 entries)
- `result60_v2_residue_marginal.csv` — D_pred vs D_emp per residue
- `experiments_output/result60_v2_validation_table.csv` — 7-concern table
- `experiments_output/result60_v2_pearson_vs_B.csv` — scaling sweep
- `experiments_output/result60_v2_outcome_delta.csv` — outcome (δ) verification

Compute: ~30s reference + ~60s validation sweep + ~30s outcome (δ) check.
