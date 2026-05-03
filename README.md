# Collatz residue-class structural analysis

**Status (2026-05-02):** Bridge result to Tao 2022 documented. Three findings consolidated: (a) prefix-decomposition theorem at modular resolutions k = 4..14; (b) `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))` with slope = 1.000 ± 0.005 at K_h = 3/log(4/3) across two independent observables (σ and first-passage), four modular resolutions, and two data scales; (c) qx+1 Cramér convergence law at q ∈ {5, 7, 9, 11} with q=5 match to 0.01%.

---

## TL;DR — what's here

For odd integers n, the total stopping time σ(n) under the Collatz map has structure indexed by the residue r = n mod 2^k via a *deterministic prefix*. Tracking the symbolic state (a, c) such that orbit value = a·m + c (where m is the integer tail of n) under the Collatz iteration until a becomes odd terminates at a_final ∈ {3^j : 1 ≤ j ≤ k}. The 2^(k−1) odd residue classes therefore collapse onto exactly k distinct conditional distributions of σ.

α_det(r) := prefix_steps(r) + K_h · log(a_final(r) / 2^k), where K_h = 3/log(4/3) ≈ 10.4282, predicts:

1. **σ-intercepts per class** (the original definition; verified at k = 4..12).
2. **First-passage-time means per class** for the orbit reaching ≤ f(N) for *any* threshold f, with offset matching Tao 2022's (5.15) leading term `K_h · log(N/f(N))` to ≤ 1 step (1%-trimmed mean) across 40 verification cells.

The bridge to Tao 2022 is structurally `s_mean(r; f) ≈ α_det(r) + K_h · log(N/f(N)) + ε`, where ε is small, observable-dependent, and stable across modular resolution and N.

---

## How to navigate this repo

### If you want the headline result
1. Read [`writeup.md`](writeup.md) — Result 1, Result 3, and the subsection "α_det predicts mean first-passage time and matches Tao (5.15) at the per-class level" within Result 3.
2. Then [`tao_bridge_findings.md`](tao_bridge_findings.md) for the TA.1/TA.2/TA.3 tightening of ε(N).

### If you want to verify a specific claim
- Each substantive claim in `writeup.md` traces to a `findings.md` entry and one or more numbered experiments. Use the experiment index below to locate the script. Run it; CSV outputs land in `experiments_output/`.

### If you want to check the prior-art positioning
- [`literature_check.md`](literature_check.md) — audit identifying Terras 1976 Lemma 4 as the asymptotic predecessor of the prefix decomposition, plus connections to Sinai 2003, Tao 2022, Bonacorsi & Bordoni 2026.

### If you want to see the audit trail with sanity-check protocol applied per finding
- [`findings.md`](findings.md) — chronological log, append-only.
- [`agent2_findings.md`](agent2_findings.md) — trajectory-measure deep dive (v=4/v=10 spike mechanism, q=5 trajectory v, MGF preservation across q).
- [`compute_threads_findings.md`](compute_threads_findings.md) — σ-record extension to OEIS A006877 b-file, prefix-tail mechanism analysis, q=5 cycle search.

### If you want the qx+1 generalization (companion to the 3x+1 work)
- The qx+1 Cramér convergence-rate result lives separately. Code is mixed in among the experiments below (10–22 range and 28/29 range). For prior consolidations see `findings.md` entries dated 2026-05-01 onward and the auto-memory file `project_collatz_qx1.md` (in the user's external memory).

---

## Documents

| File | Contents |
|---|---|
| [`writeup.md`](writeup.md) | Canonical result document. Result 1 (slope universality + non-monotone β oscillation), Result 2 (tail shape), Result 3 (prefix decomposition + Tao bridge subsection), Related Work (B&B comparison + Pathfinder caveat), Limitations. |
| [`findings.md`](findings.md) | Append-only chronological audit trail. Every empirical finding gets sanity-check protocol entries (sampling bias / definition / finite-N / parity / numerical precision). ~600 lines. |
| [`agent2_findings.md`](agent2_findings.md) | Trajectory-measure characterization: q=3 trajectory v moments, MGF preservation across q, m mod 32/2048/131072 pushforward (mechanism for v=4/v=10 spikes), q=5 unconditional trajectory v. |
| [`compute_threads_findings.md`](compute_threads_findings.md) | σ-record class-fraction analysis (T1.1/T1.2/T1.5/T1.6/TB.2). Gaussian-tail Gumbel mechanism for prefix-class σ-record fractions, replacing earlier exponential-θ guess. |
| [`tao_bridge_findings.md`](tao_bridge_findings.md) | TA.1 N-stability of σ structural offset (constant ≈ −2.45 across N = 2²⁵..2³²), TA.2 trim-quantile sweep (q* = 1.18% drives gap to 0 at √N), TA.3 parametric fit (gap ≈ −2.35 + 0.486·log(threshold)). |
| [`closed_form_findings.md`](closed_form_findings.md) | Closed-form derivations for the bridge structural constants. ⟨α_det⟩ = log(6)/log(4/3) DERIVED exactly. ε(σ) and slope-on-log(threshold) ruled out as having clean closed forms; trace back to either Lagarias trajectory measure (open) or finite-N μ_β characterization (TA.1 follow-up). |
| [`literature_check.md`](literature_check.md) | Prior-art audit. Terras 1976 Lemma 4, Sinai 2003, Lagarias 1985, Tao 2022, B&B 2026. |
| [`one_sheet_lin.py`](one_sheet_lin.py) / [`one_sheet_yosef.py`](one_sheet_yosef.py) | PDF generators for one-sheet summaries (Lin: 3x+1; Yosef: qx+1). |

---

## Experiment index (organized by theme)

Numbered by creation order; some numbers collide because the project ran two parallel agents. Filename disambiguates.

### Stage 1–4 pipeline (original Bayesian fit)
| Script | Purpose |
|---|---|
| [`generate.py`](generate.py) | Numba memoized σ / syracuse / max_excursion / residues for n ∈ [1, N]. Outputs `data/main_N{N}.parquet`. |
| [`generate_q.py`](generate_q.py) | qx+1 generalization: writes `data/q_main_q{q}_N{N}.parquet`. |
| [`analyze.py`](analyze.py) | Stage 2 EDA: σ vs log n by mod-16 class, v-distribution, residual tails. |
| [`stage3_prep.py`](stage3_prep.py) | Stage 3 input: odd-only filter, class index, uniform stratified subsample. |
| [`fit.py`](fit.py) | Stage 3 hierarchical Stan fit (k=6, k=10). Outputs to `fits/{tag}/`. |
| [`diagnose.py`](diagnose.py) | Stage 4 posterior summary, GPD on tails, posterior tail probabilities. Outputs to `stage4_results/{tag}/`. |

### Per-class structure (3x+1)
| # | Script | Purpose |
|---|---|---|
| 01 | [`experiments/01_alpha_decomposition.py`](experiments/01_alpha_decomposition.py) | Per-class OLS α(r) vs predicted α_det(r) at given k. |
| 02 | [`experiments/02_moment_universality.py`](experiments/02_moment_universality.py) | Higher per-class moments (variance, skew, kurtosis) vs prefix prediction. |
| 03 | [`experiments/03_n_scaling.py`](experiments/03_n_scaling.py) | μ_β scaling N ∈ {2²⁰..2²⁵}. |
| 05 | [`experiments/05_cfinal_ks_analysis.py`](experiments/05_cfinal_ks_analysis.py) | Within-a_final c_final substructure via KS tests. |
| 07 | [`experiments/07_anderson_darling.py`](experiments/07_anderson_darling.py) | Distributional clustering of per-class σ residuals via Anderson-Darling. |
| 08 | [`experiments/08_all_n_decomposition.py`](experiments/08_all_n_decomposition.py) | Decomposition extended to all n (odd ∪ even). |
| 09 | [`experiments/09_multi_stat_decomposition.py`](experiments/09_multi_stat_decomposition.py) | Decomposition for σ, syracuse, odd_steps, even_steps, log(max_excursion). |
| 24 | [`experiments/24_k_sweep_alpha_decomposition.py`](experiments/24_k_sweep_alpha_decomposition.py) | k-sweep at N=2²⁷ for k ∈ {4..12}; noise-floor ratio band. |

### B&B NB GLM replication
| # | Script | Purpose |
|---|---|---|
| 04 | [`experiments/04_head_to_head_nb_glm.py`](experiments/04_head_to_head_nb_glm.py) | Frequentist NB GLM head-to-head (M0..M4). |
| 06 | [`experiments/06_bb_replication.py`](experiments/06_bb_replication.py) | Bayesian NB GLM via cmdstanpy NUTS. |
| 06b | [`experiments/06b_bb_pathfinder.py`](experiments/06b_bb_pathfinder.py) | Pathfinder VI fallback (used when NUTS multi-chain locked at N_train=500K). |
| — | [`experiments/nb2_glm.stan`](experiments/nb2_glm.stan) | Shared Stan model. |

### Trajectory measure
| # | Script | Purpose |
|---|---|---|
| 15 / 15b | [`experiments/15_step_variance.py`](experiments/15_step_variance.py), [`experiments/15b_step_variance_unconditional.py`](experiments/15b_step_variance_unconditional.py) | Conditional vs unconditional v variance. |
| 25 | [`experiments/25_trajectory_measure.py`](experiments/25_trajectory_measure.py) | High-resolution v-distribution at N_start=10⁸, T=200; v=4 and v=10 spike characterization. |
| 27 | [`experiments/27_m_residue_pushforward.py`](experiments/27_m_residue_pushforward.py) | m mod 32 / 2048 / 131072 pushforward (mechanism for v=4/v=10/v=16 spikes). Agent 2. |
| 28 | [`experiments/28_per_octave_trajectory_E_v.py`](experiments/28_per_octave_trajectory_E_v.py) | Per-octave trajectory E[v] for the K(E[v]) closed-form prediction of β_local. |
| 28 | [`experiments/28_q5_trajectory_measure.py`](experiments/28_q5_trajectory_measure.py) | q=5 trajectory v-distribution on convergent orbits. Agent 2. |
| 29 | [`experiments/29_v_step_correlation.py`](experiments/29_v_step_correlation.py) | Lag-1 autocorrelation of v along Syracuse trajectories. |

### β oscillation and N-extension
| # | Script | Purpose |
|---|---|---|
| 26 | [`experiments/26_mu_beta_n_extension.py`](experiments/26_mu_beta_n_extension.py) | Streaming OLS at N up to 2³². Non-monotone β oscillation finding. |
| 27 | [`experiments/27_beta_oscillation_diagnostic.py`](experiments/27_beta_oscillation_diagnostic.py) | Per-octave β_local + top-K outlier exclusion (record-σ hypothesis test). |

### First-passage / Tao bridge
| # | Script | Purpose |
|---|---|---|
| 23 | [`experiments/23_sigma_fiber_cardinality.py`](experiments/23_sigma_fiber_cardinality.py) | σ-fiber cardinality (Avenue A diagnostic; cryptographic hardness ruled out). |
| 30 | [`experiments/30_first_passage_alpha_det.py`](experiments/30_first_passage_alpha_det.py) | First Spearman ρ = 1.0 finding for s_median vs α_det at k=8. |
| 31 | [`experiments/31_first_passage_replication.py`](experiments/31_first_passage_replication.py) | Replication at k=8/10/12 × 4 thresholds. |
| 32 | [`experiments/32_alpha_det_K_calibration.py`](experiments/32_alpha_det_K_calibration.py) | K-recalibration test for s_median (slope < 1 mechanism diagnosis). |
| 33 | [`experiments/33_alpha_det_K_calibration_mean.py`](experiments/33_alpha_det_K_calibration_mean.py) | s_mean version: slope = 1 at K_h with raw mean; trim-1% match to Tao. |
| 34 | [`experiments/34_alpha_det_K_calibration_mean_k_sweep.py`](experiments/34_alpha_det_K_calibration_mean_k_sweep.py) | k=8/10/12 mean replication. |
| 35 | [`experiments/35_alpha_det_full_bridge.py`](experiments/35_alpha_det_full_bridge.py) | Full bridge: 4 k × 5 observables (σ + 4 first-passage thresholds) × 2 N. 40-cell verification. |
| 36 | [`experiments/36_TA1_sigma_offset_N_sweep.py`](experiments/36_TA1_sigma_offset_N_sweep.py) | TA.1: σ offset N-stability at N = 2²⁵..2³². |
| 37 | [`experiments/37_TA2_trim_quantile_sweep.py`](experiments/37_TA2_trim_quantile_sweep.py) | TA.2: trim quantile sweep finding q* = 1.18%. |
| 37 | [`experiments/37_alpha_det_k16_verification.py`](experiments/37_alpha_det_k16_verification.py) | Higher-k verification at k=16. Agent 2. |
| 38 | [`experiments/38_TA3_parametric_fit.py`](experiments/38_TA3_parametric_fit.py) | TA.3: parametric fit `gap ≈ −2.35 + 0.486 · log(threshold)`. |
| 39 | [`experiments/39_overshoot_at_first_passage.py`](experiments/39_overshoot_at_first_passage.py) | First-passage overshoot ⟨log(f/v*)⟩ — empirical ≈ 0.298 nats, constant in f. Rules out the K_h·⟨overshoot⟩ explanation for the slope on log(f). |
| 40 | [`experiments/40_K_eff_decomposition.py`](experiments/40_K_eff_decomposition.py) | Parity decomposition of K_eff into odd-σ slope + halving compensation; user's hypothesis K_eff = 9.31 + 0.63 ruled out empirically. |

### qx+1 generalization (companion)
| # | Script | Purpose |
|---|---|---|
| 10 | [`experiments/10_q_decomposition.py`](experiments/10_q_decomposition.py) | qx+1 prefix decomposition at k=6 for q ∈ {3, 5, 7, 11}. |
| 10b | [`experiments/10b_q_partial_correlation.py`](experiments/10b_q_partial_correlation.py) | Partial correlation diagnostics (j-slope vs log(m)). |
| 12 | [`experiments/12_q_convrate_analytical.py`](experiments/12_q_convrate_analytical.py) | log(conv_rate) vs j slope fits. |
| 13 | [`experiments/13_cross_q_unification.py`](experiments/13_cross_q_unification.py) | C ≈ 5/2 universal-multiplier hypothesis test (rejected at q=11). |
| 14 | [`experiments/14_conv_rate_vs_N.py`](experiments/14_conv_rate_vs_N.py) | conv_rate(N) decay exponent. |
| 16 | [`experiments/16_cramer_root.py`](experiments/16_cramer_root.py) | Exact Cramér root: `q^(−θ) = 2^(1−θ) − 1`. |
| 17 | [`experiments/17_cramer_dual_verification.py`](experiments/17_cramer_dual_verification.py) | Dual j-slope and N-decay verification of θ(q). |
| 18 | [`experiments/18_q7_x_binning_diagnostic.py`](experiments/18_q7_x_binning_diagnostic.py) | Pooled-X vs per-class diagnostic at q=7. |
| 19 | [`experiments/19_bahadur_rao.py`](experiments/19_bahadur_rao.py) | Bahadur-Rao 1/√L sub-exponential prefactor test (rejected). |
| 20 | [`experiments/20_m_selection_test.py`](experiments/20_m_selection_test.py) | m-selection partial correlation diagnostic. |
| 21 | [`experiments/21_two_term_fit.py`](experiments/21_two_term_fit.py) | Two-term fit `f(X) = A·X^(−θ) + B` (B → 0; rejected). |
| 22 | [`experiments/22_q5_cycle_detection.py`](experiments/22_q5_cycle_detection.py) | Floyd cycle detection at q=5; non-trivial cycle landings 0.12%. |
| 29 | [`experiments/29_qx1_cycle_classification.py`](experiments/29_qx1_cycle_classification.py) | qx+1 cycle classification at q ∈ {5, 7, 11, 13}. Agent 2. |
| 36 | [`experiments/36_q5_fourth_cycle_search.py`](experiments/36_q5_fourth_cycle_search.py) | q=5 fourth-cycle search. Agent 2. |

### Auxiliary
| # | Script | Purpose |
|---|---|---|
| 30 | [`experiments/30_sigma_records_prefix_analysis.py`](experiments/30_sigma_records_prefix_analysis.py) | σ-record class-fraction Gaussian-Gumbel analysis. Agent 2. |

---

## Data files

`data/main_N{N}.parquet`:

| N (= 2^k or scientific) | rows | size |
|---|---|---|
| 2²⁰ = 1,048,576 | 1.0M | 7.5 MB |
| 4,194,304 | 4.2M | 27 MB |
| 8,388,608 | 8.4M | 52 MB |
| 10,000,000 | 10M | 62 MB |
| 16,777,216 | 16.8M | 102 MB |
| 33,554,432 = 2²⁵ | 33.6M | 201 MB |
| 134,217,728 = 2²⁷ | 134M | 786 MB |

Schema: `n, sigma, syracuse, odd_steps, even_steps, max_excursion, is_record, res_mod_16, res_mod_64, res_mod_256`.

`data/q_main_q{q}_N{N}.parquet`: qx+1 versions for q ∈ {3, 5, 7, 9, 11, 13} at varying N up to 10⁹.

Larger σ caches at N ∈ {2²⁸, 2³⁰, 2³²} are built in-memory by experiment 36 and not persisted (~1, 4, 17 GB int32 respectively).

---

## Reproduction smoke check

If these three pass, the work is intact:

1. **Data generation:** `python generate.py --N 1048576` should finish in ~1.5 s. σ at n=27 should be 111.
2. **Prefix algorithm:** Run by hand on residue r = 21 starting from state (a=64, c=21). Expected: 7 steps, terminating at (a_final=3, c_final=1).
3. **Bridge result:** `python experiments/35_alpha_det_full_bridge.py` should report slope at K_h ≈ 1.000 ± 0.005 across all 40 cells, with offset gaps matching the table in `tao_bridge_findings.md`.

For the headline N-extension finding: `python experiments/36_TA1_sigma_offset_N_sweep.py` reports gap ≈ −2.45 across N = 2²⁵ → 2³², stable to 0.01.

Compute requirements:
- All experiments at N ≤ 2²⁷ run in seconds-to-minutes on a 16-thread CPU.
- N=2³⁰ sigma cache: ~14 s, 4 GB RAM.
- N=2³² sigma cache: ~55 s, 17 GB RAM.

---

## Open follow-ups

Not load-bearing for the bridge claim, but adjacent and worth pursuing:

- **Closed-form derivation of the −2.35 structural constant** (TA.3 intercept). The N-stability data (TA.1) supports it being a structural invariant of the σ distribution; analytical form unresolved.
- **Closed-form derivation of the +0.486 ≈ 1/2 slope on log(threshold)** (TA.3). Hints at a √f scaling but the underlying mechanism is not isolated.
- **Trim-quantile interpretation at larger N.** TA.2's q* = 1.18% does not match a clean log^(−c) N exponent; might be an N=2²⁷ artifact. Re-test at N=2³².
- **q=5 / q=7 first-passage analog.** Does the bridge `s_mean ≈ α_det^(q) + K_q · Δlog` hold for qx+1 with K_q derived from θ(q)? Generalizes the bridge across the qx+1 family.
- **Bonacorsi-side HMC validation** at full N=10⁷ on the NB GLM replication (deferred to N. Bonacorsi at Columbia).

---

## Outreach packages on Desktop

- `collatz_for_lin_2026-05-01/` and `.zip` — Lin (Maryland) review package: writeup, findings, one-sheet PDF, data sample.
- `collatz_for_bonacorsi_2026-05-01_v4_pathfinder/` — B&B framework replication artifacts.
- `collatz_qx1_for_yosef_2026-05-01/` — qx+1 Cramér derivation package for the probability-theory audience.

These predate the Tao bridge result. Re-zipping with the updated `writeup.md` and `tao_bridge_findings.md` would refresh them.
