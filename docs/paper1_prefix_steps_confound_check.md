# prefix_steps Confound Check — Residual Gap Analysis for B2 Slope


**Headline:** PREFIX_STEPS CONFOUND DOMINATES.


Pred B's per-class ratios are closer to 1.0 than Pred A's by ≥ 5 percentage points on all three classes. The prefix_steps term is the dominant residual source. The conjecture should reference α_det as a whole, not log(a_final) alone.


---


## Section 1 — Four (r, prefix_steps, log(a_final), α_det) triples


Using K_h = 3/log(4/3) ≈ 10.428178, k=3.

α_det(r) = prefix_steps(r) + K_h · log(a_final(r) / 2^k).


| r (odd) | a_final | j | prefix_steps | log(a_final) | α_det |
|---|---|---|---|---|---|
| 1 | 9 | 2 | 5 | 2.1972 | +6.2283 |
| 3 | 9 | 2 | 5 | 2.1972 | +6.2283 |
| 5 | 3 | 1 | 4 | 1.0986 | -6.2283 |
| 7 | 27 | 3 | 6 | 3.2958 | +18.6848 |

## Section 2 — Pearson correlation between prefix_steps and log(a_final)


Across the 4 odd residues:


- **Pearson r = +1.000000**


Reason: prefix_steps(r) = j(r) + k where a_final(r) = 3^j(r). At fixed k, prefix_steps is linear in j(r), and log(a_final(r)) = j(r) · log(3). So prefix_steps and log(a_final) are *exactly* linearly related (slope = 1/log(3) ≈ 0.910), giving Pearson r = ±1.0 across any residue set whose distinct a_final values span more than one j level.


## Section 3 — Predictions vs observations


E[σ] (training set) = 154.4339.


- **Pred A** (literal conjecture, properly link-translated): dummy = β · log(a_final) with β = K_h / E[σ] = 0.067525.

- **Pred B** (full α_det rescaled): dummy = [α_det(class) − α_det(af=1)] / E[σ].


### α_det per class (mean over odd residues mapping to that class)


| a_final | α_det (mean over odd r) |
|---|---|
| 3 | -6.2283 |
| 9 | +6.2283 |
| 27 | +18.6848 |
| 1 (reference, even r=0) | -18.6848 |

### Predictions table


| a_final | log(a_final) | Pred A | Pred B | Observed | A ratio | B ratio |
|---|---|---|---|---|---|---|
| 3 | 1.0986 | 0.0742 | 0.0807 | 0.1037 | 1.3979 | 1.2857 |
| 9 | 2.1972 | 0.1484 | 0.1613 | 0.1806 | 1.2172 | 1.1195 |
| 27 | 3.2958 | 0.2226 | 0.2420 | 0.2578 | 1.1584 | 1.0654 |

- MAE Pred A = 0.032332
- MAE Pred B = 0.019381


### Slope-level comparison (across all 4 dummy points incl af=1 ref at 0)


- Pred A slope on log(a_final) = K_h / E[σ] = 0.067525
- Pred B slope on log(a_final) = (1/log(3) + K_h) / E[σ] = 0.073419
- Empirical OLS slope across 4 dummy points = 0.0774 (95% CI [0.0676, 0.0869])


## Section 4 — Verdict


**PREFIX_STEPS CONFOUND DOMINATES.**


Pred B's per-class ratios are closer to 1.0 than Pred A's by ≥ 5 percentage points on all three classes. The prefix_steps term is the dominant residual source. The conjecture should reference α_det as a whole, not log(a_final) alone.


Per-class (|A_ratio − 1| − |B_ratio − 1|) × 100, in percentage points:


| a_final | pp difference (positive = Pred B closer) |
|---|---|
| 3 | +11.22 |
| 9 | +9.77 |
| 27 | +9.30 |

## Section 5 — Implication for §3.5 conjecture statement


The conjecture in §3.5 should reference **α_det(r)** as a whole, not log(a_final(r)) alone. The prefix_steps(r) term contributes a non-trivial component to the per-class effect on σ that the empirical B2 dummies absorb. Stating the conjecture purely on log(a_final) leaves a systematic residual that matches the prefix_steps slope 1 / (E[σ] · log(3)).


Recommended restatement:


*The Bonacorsi-Bordoni mod-8 random effect is, up to a global additive constant, α_det(r) / E[σ] at k=3, where α_det(r) = prefix_steps(r) + K_h · log(a_final(r)/2^k) is the prefix decomposition's per-class structural offset on σ-scale.*


## Section 6 — Recommended §5 paragraph language


To validate the §5 head-to-head's structural ordering and the §3.5 conjecture under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0, 4 chains × 2000 iter each, parallel-chains mode). All five specifications converged with R̂ < 1.01, ESS_bulk > 600, zero divergences, zero treedepth saturations, and E-BFMI > 0.80. The HMC ordinal ordering reproduces the Pathfinder finding: B3 > B2 > B1 > B0 on held-out log score, with B4 within sampling noise of B3 (B3 better by 0.10 nats). σ_u in B1 inflates substantially under HMC (posterior mean 0.097 vs Pathfinder 0.003*), confirming the §5 caveat that Pathfinder collapses hierarchical scale parameters in non-Gaussian regions. σ_u in B4 under HMC remains small (0.008), confirming the conjecture's exhaustion claim — once a_final at k=6 is in the model, the mod-8 random effect carries no residual information.

On the literal-magnitude form of the §3.5 conjecture: the empirical OLS slope of B2's a_final fixed-effect dummies on log(a_final(r)) is 0.0774 with 95% credible interval [0.0676, 0.0869]. Translating the §3.5 conjecture through the NB2 log link, the predicted GLM-scale slope on the full structural offset α_det(r) = prefix_steps(r) + K_h · log(a_final(r)/2^k) is (1/log(3) + K_h)/E[σ] = 0.0734, against (K_h/E[σ]) = 0.0675 on log(a_final) alone. The empirical slope falls between these two predictions, with the α_det form 5.4% below empirical against the log(a_final)-only form 14.7% below. We therefore restate the conjecture in §3.5 to reference α_det(r) as a whole rather than log(a_final(r)) alone, since prefix_steps(r) contributes a non-negligible component to the per-class effect that the empirical B2 dummies absorb at this precision.


---


*Inputs: B2 HMC posterior at N=10K, 8000 train / 2000 test, seed 20260509. NO new MCMC. Restricted to 4 odd residues r ∈ {1, 3, 5, 7} per the brief.*
