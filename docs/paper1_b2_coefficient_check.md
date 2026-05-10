# B2 Coefficient and log(n) Correlation Check — Pre-Submission Verdict

**Headline:** CONJECTURE FALSIFIED IN MAGNITUDE.


B2's a_final dummy coefficients exclude a slope of 1.0 on log(a_final) and are consistent with the small B1 magnitude. The literal-magnitude form of the conjecture is rejected; the directional/monotone form survives.

---


## Section 1 — Test 1: log(n) ↔ log(a_final(r mod 8)) correlation


- Pearson r = +0.004145 (p = 7.1090e-01), 95% CI [-0.0178, +0.0261]
- Spearman ρ = +0.001216 (p = 9.1337e-01)

**Verdict: INDEPENDENT.** |r|<0.02; explanation (a) parameterization artifact is RULED OUT. log(n) and log(a_final(r mod 8)) are essentially independent.


## Section 2 — Test 2: B2 a_final dummy coefficients


B2's design uses categorical dummies on a_final (reference = af=1, j=0). Posterior coefficients on the 3 non-reference levels:


| level | a_final | j | log(a_final) | β posterior mean | 95% CI | conjecture predicts |
|---|---|---|---|---|---|---|
| af = 3 | 3^1 | 1 | 1.0986 | **+0.1037** | [+0.0762, +0.1313] | +1.0986 (✗ excludes 1.0×log(af)) |
| af = 9 | 3^2 | 2 | 2.1972 | **+0.1806** | [+0.1534, +0.2077] | +2.1972 (✗ excludes 1.0×log(af)) |
| af = 27 | 3^3 | 3 | 3.2958 | **+0.2578** | [+0.2234, +0.2919] | +3.2958 (✗ excludes 1.0×log(af)) |

Implied per-draw OLS slope on log(a_final) (4 dummy values, a_final=1 reference fixed at 0):

- **Posterior mean slope = +0.0774**, 95% CI [+0.0676, +0.0869]
- Conjecture's literal magnitude predicts slope = 1.0
- B1's smoke-test "falsified" form predicts slope ≈ 0.085

**Verdict: LITERAL-MAGNITUDE FALSIFIED.** 95% CI [+0.0676, +0.0869] contains 0.085 and EXCLUDES 1.0. The literal conjecture's magnitude is rejected; the smoke-test verdict stands.


## Section 3 — Test 3: β_log_n in B0 vs B2


| spec | β_log_n posterior mean | SD | 95% CI |
|---|---|---|---|
| B0 (log(n) only) | 0.072791 | 0.004435 | [0.063933, 0.081187] |
| B2 (+ a_final dummies) | 0.073357 | 0.004522 | [0.064371, 0.082093] |

- Difference (B2 − B0) = +0.000566, pooled SD = 0.006333, |z| = 0.089

**Verdict: ESSENTIALLY UNCHANGED.** |z| = 0.09 ≤ 2: B2's β_log_n is essentially identical to B0's. log(n) and a_final dummies are nearly orthogonal in the training data — no signal redistribution.


## Section 4 — Test 4: Parameterization-corrected u_r prediction


For each residue mod 8, compute the predicted random effect under the literal conjecture, corrected for any colinearity with log(n):

```
u_r_predicted = log(a_final(r)) - β_log_n × (mean(log(n) | r) - global mean(log(n)))
```

- Global mean(log(n)) = 15.1235

| r | mean(log(n)\|r) | log(a_final) | u_r predicted (corrected) | u_r actual (B1) |
|---|---|---|---|---|
| 0 | 15.1303 | 0.0000 | -0.0005 | -0.1337 |
| 1 | 15.1164 | 2.1972 | +2.1977 | +0.0427 |
| 2 | 15.1421 | 1.0986 | +1.0973 | -0.0352 |
| 3 | 15.0814 | 2.1972 | +2.2003 | +0.0540 |
| 4 | 15.0921 | 1.0986 | +1.1009 | -0.0454 |
| 5 | 15.1330 | 1.0986 | +1.0979 | -0.0142 |
| 6 | 15.1236 | 2.1972 | +2.1972 | +0.0350 |
| 7 | 15.1660 | 3.2958 | +3.2927 | +0.1188 |

- Range (predicted, parameterization-corrected) = 3.2932
- Range (actual B1) = 0.2525
- Pearson r (demeaned predicted vs actual) = +0.9873 (p = 5.1014e-06)
- Spearman ρ = +0.9286 (p = 8.6297e-04)

If the parameterization correction recovers the literal magnitude, the predicted range should be ≈ 3.30 (= log(27)) and the demeaned correlation should be high. If the correction has no effect, the predicted range stays at log(27) (since the colinearity is small) and the actual range stays at 0.25 — meaning the magnitude gap is real, not a parameterization artifact.


## Section 5 — Combined Verdict


**CONJECTURE FALSIFIED IN MAGNITUDE.**


B2's a_final dummy coefficients exclude a slope of 1.0 on log(a_final) and are consistent with the small B1 magnitude. The literal-magnitude form of the conjecture is rejected; the directional/monotone form survives.


## Section 6 — Recommended §5 Paragraph (draft for editorial review)


To validate the §5 head-to-head's structural ordering and the §3.5 conjecture under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0, 4 chains × 2000 iter each, parallel-chains mode). All five specifications converged with R̂ < 1.01, ESS_bulk > 600, zero divergences, zero treedepth saturations, and E-BFMI > 0.80. The HMC ordinal ordering reproduces the Pathfinder finding: B3 > B2 > B1 > B0 on held-out log score, with B4 within sampling noise of B3 (B3 better by 0.10 nats). σ_u in B1 inflates substantially under HMC (posterior mean ≈ 0.097 vs Pathfinder 0.003*), confirming the §5 caveat that Pathfinder collapses hierarchical scale parameters. σ_u in B4 under HMC remains small (0.008), confirming that the mod-8 random effect carries no residual information once a_final at k=6 is in the model.

On the literal-magnitude form of the §3.5 conjecture (u_r = log(a_final(r)) + const at k=3): B2's posterior coefficients on the a_final fixed effects, when regressed against log(a_final), give a slope of +0.077 with 95% credible interval [+0.068, +0.087], which excludes 1.0. The literal magnitude is rejected at this scale. The DIRECTIONAL form survives — u_r is monotone-increasing in log(a_final(r)) — but the slope is materially smaller than log(3). We therefore restate the conjecture as a directional claim: *the mod-8 random effect is monotone in the prefix-decomposition's log(a_final(r)) at k=3, with a slope to be determined empirically.* The structural origin claim — that the random effect's signal traces to the prefix's terminal a_final — survives; the algebraic identity at unit slope does not.

---


*Inputs: B0/B1/B2 HMC posteriors at N=10K, 8000 train / 2000 test, seed 20260509. NO new MCMC. Tests 1–4 on existing chains.*
